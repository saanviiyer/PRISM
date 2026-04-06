import os
import csv
import random
import torch
import esm
import numpy as np
from Bio.Seq import Seq
from datetime import datetime
from scipy.spatial.distance import cosine

# ---- Mock Reverse Translation ----
CODON_MAP = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGC',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
    '*': 'TAA'
}

def mock_reverse_translate(aa_seq):
    clean_seq = aa_seq.upper().replace('\n', '')
    return "".join([CODON_MAP.get(aa, 'NNN') for aa in clean_seq if aa in CODON_MAP])

# ---- ESM Setup ----
def load_esm_model():
    print("Loading ESM-2 model (this may take a moment)...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    return model, alphabet, batch_converter

def get_esm_embedding(aa_seq, model, alphabet, batch_converter):
    aa_seq = aa_seq[:1022]
    data = [("seq", aa_seq.upper())]
    _, _, tokens = batch_converter(data)
    with torch.no_grad():
        results = model(tokens, repr_layers=[33], return_contacts=False)
    token_representations = results["representations"][33]
    return token_representations[0, 1 : len(aa_seq) + 1].mean(0).numpy()

# ---- Dynamic Profile Discovery ----
def discover_profiles(base_dir):
    profiles_to_test = {}
    if not os.path.exists(base_dir):
        print(f"Warning: Directory '{base_dir}' not found.")
        return profiles_to_test
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.txt') and not file.startswith('.'):
                category = os.path.basename(root)
                file_name_clean = file.replace('.txt', '')
                signature_label = f"{category} ({file_name_clean})"
                profiles_to_test[signature_label] = os.path.join(root, file)
    return profiles_to_test

# ---- Mutation Logic ----

def detect_signature_type(profile_path):
    """
    Detects whether a profile file is SBS, DBS, or ID format by inspecting
    the first data row's mutation context key.

    SBS format:  A[C>A]A   — trinucleotide context with brackets
    DBS format:  AC>CA     — two-base change with '>' and no brackets
    ID format:   1:Del:C:0 — colon-separated indel descriptor
    """
    with open(profile_path, 'r') as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split()
            if parts:
                key = parts[0]
                if '[' in key and ']' in key:
                    return 'SBS'
                elif ':' in key:
                    return 'ID'
                elif '>' in key:
                    return 'DBS'
    return 'SBS'  # fallback

def load_mutation_profile(profile_path, column_name=None):
    """
    Loads a mutation profile from a .txt file, auto-detecting SBS / DBS / ID format.

    For SBS and DBS: expects multiple genome columns in the header.
      column_name selects which genome version to use (e.g. 'GRCh38').
      Partial matching is supported: 'GRCh38' matches 'SBS8_GRCh38'.

    For ID: typically has only one data column. column_name still attempted
      but falls back gracefully if not found.

    Returns a dict: { context_key: probability, ... }
    Also returns the detected signature type string: 'SBS', 'DBS', or 'ID'.
    """
    sig_type = detect_signature_type(profile_path)
    sigs = {}

    with open(profile_path, 'r') as f:
        header = f.readline().strip().split()

        # Resolve which column index to read probabilities from
        col_idx = 1  # default: first data column
        if column_name is not None:
            if column_name in header:
                col_idx = header.index(column_name) + 1
            else:
                matches = [h for h in header if column_name in h]
                if matches:
                    col_idx = header.index(matches[0]) + 1
                elif len(header) >= 1:
                    # ID files often have only 1 header col — use col 1 silently
                    if sig_type != 'ID':
                        print(f"    WARNING: Column '{column_name}' not found in "
                              f"{os.path.basename(profile_path)}. "
                              f"Available: {header}. Using column 1.")

        for line in f:
            parts = line.strip().split()
            if len(parts) > col_idx:
                try:
                    sigs[parts[0]] = float(parts[col_idx])
                except ValueError:
                    continue

    return sigs, sig_type


# ---- Mutation Applicators (one per signature type) ----

def apply_sbs_mutation(nt_seq, profile):
    """
    SBS — Single Base Substitution.
    Context key format: A[C>A]A
      - Characters 0,2,6 = left_base, ref_base, right_base  → trinucleotide context
      - Character 4      = alt_base (the substitution)
    Scans every position in the sequence for matching trinucleotide contexts,
    then samples one mutation weighted by the signature probabilities.
    """
    seq_list = list(nt_seq.upper())
    possible_mutations = []
    for idx in range(len(seq_list) - 2):
        context = "".join(seq_list[idx:idx + 3])
        for mut_sig, prob in profile.items():
            if len(mut_sig) < 7 or prob <= 0:
                continue
            orig_tn = mut_sig[0] + mut_sig[2] + mut_sig[6]
            if context == orig_tn:
                possible_mutations.append({
                    'index': idx + 1,
                    'new_base': mut_sig[4],
                    'weight': prob
                })
    if not possible_mutations:
        return nt_seq, False
    weights = [m['weight'] for m in possible_mutations]
    total = sum(weights)
    probs = [w / total for w in weights]
    chosen = np.random.choice(possible_mutations, p=probs)
    seq_list[chosen['index']] = chosen['new_base']
    return "".join(seq_list), True


def apply_dbs_mutation(nt_seq, profile):
    """
    DBS — Doublet Base Substitution.
    Context key format: AC>CA
      - Characters 0,1 = the two reference bases (dinucleotide)
      - Characters 3,4 = the two alternate bases
    Scans every adjacent pair of bases for a matching dinucleotide,
    then replaces both simultaneously.
    """
    seq_list = list(nt_seq.upper())
    possible_mutations = []
    for idx in range(len(seq_list) - 1):
        dinuc = seq_list[idx] + seq_list[idx + 1]
        for mut_sig, prob in profile.items():
            if len(mut_sig) < 5 or prob <= 0:
                continue
            ref_dinuc = mut_sig[0] + mut_sig[1]
            alt_dinuc = mut_sig[3] + mut_sig[4]
            if dinuc == ref_dinuc and ref_dinuc != alt_dinuc:
                possible_mutations.append({
                    'index': idx,
                    'new_bases': alt_dinuc,
                    'weight': prob
                })
    if not possible_mutations:
        return nt_seq, False
    weights = [m['weight'] for m in possible_mutations]
    total = sum(weights)
    probs = [w / total for w in weights]
    chosen = np.random.choice(possible_mutations, p=probs)
    seq_list[chosen['index']]     = chosen['new_bases'][0]
    seq_list[chosen['index'] + 1] = chosen['new_bases'][1]
    return "".join(seq_list), True


def apply_id_mutation(nt_seq, profile):
    """
    ID — Insertion/Deletion.
    Context key format: 1:Del:C:0  (size : type : base : repeat_count)
      - 'Del' removes one or more bases at a randomly chosen valid position.
      - 'Ins' inserts one or more bases at a randomly chosen position.
    Since indels shift the reading frame, small indels (1–2 bp) dominate
    real mutational signatures and are most biologically relevant here.
    Samples an indel type weighted by signature probability, then applies
    it at a random eligible position in the sequence.
    """
    seq_list = list(nt_seq.upper())
    possible_events = []
    for mut_sig, prob in profile.items():
        if prob <= 0:
            continue
        parts = mut_sig.split(':')
        if len(parts) < 2:
            continue
        size_str = parts[0]
        op = parts[1]  # 'Del' or 'Ins'
        try:
            size = int(size_str)
        except ValueError:
            size = 1
        possible_events.append({'op': op, 'size': size, 'weight': prob})

    if not possible_events or len(seq_list) < 2:
        return nt_seq, False

    weights = [e['weight'] for e in possible_events]
    total = sum(weights)
    probs = [w / total for w in weights]
    chosen = np.random.choice(possible_events, p=probs)

    op   = chosen['op']
    size = chosen['size']

    if op == 'Del':
        if len(seq_list) <= size:
            return nt_seq, False
        pos = random.randint(0, len(seq_list) - size - 1)
        del seq_list[pos:pos + size]
    else:  # Ins
        pos = random.randint(0, len(seq_list))
        insert = [random.choice('ACGT') for _ in range(size)]
        seq_list[pos:pos] = insert

    return "".join(seq_list), True


def apply_mutation_by_type(nt_seq, profile, sig_type):
    """Dispatch to the correct applicator based on detected signature type."""
    if sig_type == 'DBS':
        return apply_dbs_mutation(nt_seq, profile)
    elif sig_type == 'ID':
        return apply_id_mutation(nt_seq, profile)
    else:  # SBS (default)
        return apply_sbs_mutation(nt_seq, profile)


def apply_null_mutation(nt_seq):
    """
    NULL MODEL: Apply one mutation with uniform probability across all
    positions and all possible base substitutions. No signature bias.
    """
    seq_list = list(nt_seq.upper())
    idx = random.randint(0, len(seq_list) - 1)
    original = seq_list[idx]
    bases = [b for b in 'ACGT' if b != original]
    seq_list[idx] = random.choice(bases)
    return "".join(seq_list), True

def translate_nt_to_aa(nt_seq):
    trimmed = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)]
    return str(Seq(trimmed).translate(to_stop=True))

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        return "".join([l.strip() for l in lines if not l.startswith('>')])

# ---- Run simulation for one pair + one profile ----
def run_simulation(p1_nt_orig, p2_embedding, p1_orig_embedding,
                   profile, sig_type, model, alphabet, batch_converter,
                   max_steps, use_null=False):
    """
    Runs the mutation simulation for a single (protein pair, signature) combo.
    Dispatches to the correct mutation applicator based on sig_type.
    Returns the full trajectory of cosine distances.
    """
    initial_distance = cosine(p1_orig_embedding, p2_embedding)
    current_nt = p1_nt_orig
    trajectory = [initial_distance]

    for step in range(1, max_steps + 1):
        if use_null:
            current_nt, mutated = apply_null_mutation(current_nt)
        else:
            current_nt, mutated = apply_mutation_by_type(current_nt, profile, sig_type)

        if not mutated:
            print("  -> Simulation halted early: No valid context mutations left.")
            break

        current_aa = translate_nt_to_aa(current_nt)
        if not current_aa:
            break

        current_embedding = get_esm_embedding(current_aa, model, alphabet, batch_converter)
        dist = cosine(current_embedding, p2_embedding)
        trajectory.append(dist)

    return trajectory

# ---- CSV Output ----
def save_detailed_results_csv(all_results, output_dir):
    filepath = os.path.join(output_dir, 'prism_detailed_results.csv')
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['protein_pair', 'signature', 'sig_type', 'net_movement_toward_p2'])
        for pair_label, pair_results in all_results.items():
            for sig, (net_change, stype) in pair_results.items():
                writer.writerow([pair_label, sig, stype, f"{net_change:.6f}"])
    print(f"\n  Saved detailed results -> {filepath}")
    return filepath

def save_cross_pair_summary_csv(cross_pair_avg, sig_type_map, null_cross_avg, output_dir):
    filepath = os.path.join(output_dir, 'prism_cross_pair_summary.csv')
    sorted_cross = sorted(cross_pair_avg.items(), key=lambda x: x[1][0], reverse=True)
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank', 'signature', 'sig_type', 'avg_movement_toward_p2',
                         'std', 'n_pairs', 'beats_null_model'])
        for rank, (sig, (avg, std, n)) in enumerate(sorted_cross, 1):
            beats_null = 'YES' if avg > null_cross_avg and sig != '~~NULL (uniform random)' else 'NO'
            stype = sig_type_map.get(sig, '?')
            writer.writerow([rank, sig, stype, f"{avg:.6f}", f"{std:.6f}", n, beats_null])
    print(f"  Saved cross-pair summary -> {filepath}")
    return filepath

def save_wide_pivot_csv(all_results, output_dir):
    filepath = os.path.join(output_dir, 'prism_pivot_table.csv')
    pair_labels = list(all_results.keys())
    all_sigs = sorted(set(
        sig for pair_results in all_results.values() for sig in pair_results
    ))
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['signature'] + pair_labels + ['avg_across_pairs'])
        for sig in all_sigs:
            row_vals = [all_results[pair][sig][0] if sig in all_results[pair] else '' for pair in pair_labels]
            numeric_vals = [v for v in row_vals if v != '']
            avg = np.mean(numeric_vals) if numeric_vals else ''
            formatted = [f"{v:.6f}" if v != '' else '' for v in row_vals]
            avg_str = f"{avg:.6f}" if avg != '' else ''
            writer.writerow([sig] + formatted + [avg_str])
    print(f"  Saved pivot table      -> {filepath}")
    return filepath

# ---- Main ----
def main():
    base_path = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM'
    profile_base_directory = os.path.join(base_path, 'profiles')
    MAX_STEPS = 50

    # ----------------------------------------------------------------
    # Which genome column to use from each profile file.
    # Must match a column header in your .txt files exactly, or as a
    # substring (e.g. 'GRCh38' will match 'SBS8_GRCh38').
    # Options typically available: GRCh37, GRCh38, mm9, mm10, rn6
    # ----------------------------------------------------------------
    GENOME_COLUMN = 'GRCh38'

    # ----------------------------------------------------------------
    # Output directory — timestamped so each run is preserved separately
    # ----------------------------------------------------------------
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(base_path, f'results_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Results will be saved to: {output_dir}")
    print(f"Using genome column: '{GENOME_COLUMN}' from each profile file.")

    # ----------------------------------------------------------------
    # Define all protein pairs
    # Add or remove pairs here as needed. Format:
    #   'PAIR_LABEL': ('p1_filename.fa', 'p2_filename.fa', 'subfolder')
    # ----------------------------------------------------------------
    PROTEIN_PAIRS = {
        'PAIR1-PDK': ('A0A010RJQ9.fa', 'A0A010S4M5.fa', 'PAIR1-PDK'),
        'PAIR2-Ras': ('A0A010R5C7.fa', 'A0A010SAT4.fa', 'PAIR2-Ras'),
        'PAIR3-p53': ('A0A087VQ67.fa', 'A0A087WT22.fa', 'PAIR3-p53'),
        'PAIR4-ABC': ('A0A010Z411.fa', 'A0A011NXC8.fa', 'PAIR4-ABC'),
        'PAIR5-HR':  ('A0A009T0V7.fa', 'A0A010QZB6.fa', 'PAIR5-HR'),
    }

    # 1. Discover all mutational signature profiles
    PROFILES_TO_TEST = discover_profiles(profile_base_directory)
    if not PROFILES_TO_TEST:
        print("No profiles found. Check your 'profiles' folder path.")
        return
    print(f"Found {len(PROFILES_TO_TEST)} mutational signature profiles.")

    # 2. Load ESM-2 once (expensive — do not reload per pair)
    model, alphabet, batch_converter = load_esm_model()

    # ----------------------------------------------------------------
    # Master results table: results[pair][signature] = net_change
    # ----------------------------------------------------------------
    all_results = {}

    # 3. Loop over every protein pair
    for pair_label, (p1_file, p2_file, subfolder) in PROTEIN_PAIRS.items():
        print(f"\n{'='*70}")
        print(f"PROTEIN PAIR: {pair_label}")
        print(f"{'='*70}")

        p1_path = os.path.join(base_path, 'protein_pairs', subfolder, p1_file)
        p2_path = os.path.join(base_path, 'protein_pairs', subfolder, p2_file)

        if not os.path.exists(p1_path) or not os.path.exists(p2_path):
            print(f"  WARNING: Files not found for {pair_label}, skipping.")
            continue

        # Load and prepare sequences
        p1_aa_orig = read_fasta(p1_path)
        p2_aa      = read_fasta(p2_path)

        p1_nt_orig  = mock_reverse_translate(p1_aa_orig)
        p1_aa_sync  = translate_nt_to_aa(p1_nt_orig)

        p2_embedding       = get_esm_embedding(p2_aa,      model, alphabet, batch_converter)
        p1_orig_embedding  = get_esm_embedding(p1_aa_sync, model, alphabet, batch_converter)
        initial_distance   = cosine(p1_orig_embedding, p2_embedding)

        print(f"Initial cosine distance (P1 → P2): {initial_distance:.6f}")

        pair_results = {}

        # ---- A. Run all mutational signature profiles ----
        for signature_name, profile_path in sorted(PROFILES_TO_TEST.items()):
            profile, sig_type = load_mutation_profile(profile_path, column_name=GENOME_COLUMN)
            print(f"\n  --- Signature: {signature_name} [{sig_type}] ---")

            trajectory = run_simulation(
                p1_nt_orig, p2_embedding, p1_orig_embedding,
                profile, sig_type, model, alphabet, batch_converter,
                MAX_STEPS, use_null=False
            )

            net_change = initial_distance - trajectory[-1]
            pair_results[signature_name] = (net_change, sig_type)
            print(f"  Final distance: {trajectory[-1]:.6f} | Movement toward P2: {net_change:.6f}")

        # ---- B. Run NULL model (uniform random mutations) ----
        # Run multiple times and average to reduce stochasticity
        NULL_REPLICATES = 5
        print(f"\n  --- NULL MODEL (uniform, averaged over {NULL_REPLICATES} replicates) ---")
        null_net_changes = []
        for rep in range(NULL_REPLICATES):
            trajectory = run_simulation(
                p1_nt_orig, p2_embedding, p1_orig_embedding,
                profile=None, sig_type=None, model=model, alphabet=alphabet,
                batch_converter=batch_converter,
                max_steps=MAX_STEPS, use_null=True
            )
            null_net_changes.append(initial_distance - trajectory[-1])

        null_avg = np.mean(null_net_changes)
        null_std = np.std(null_net_changes)
        pair_results['~~NULL (uniform random)'] = (null_avg, 'NULL')
        print(f"  Null model avg movement: {null_avg:.6f} ± {null_std:.6f}")

        all_results[pair_label] = pair_results

        # ---- Per-pair leaderboard ----
        print(f"\n  LEADERBOARD — {pair_label}")
        print(f"  {'Signature':<50} {'Type':>5} | Net movement toward P2")
        print(f"  {'-'*75}")
        sorted_pair = sorted(pair_results.items(), key=lambda x: x[1][0], reverse=True)
        for rank, (sig, (change, stype)) in enumerate(sorted_pair, 1):
            better_than_null = " ✓" if change > null_avg and sig != '~~NULL (uniform random)' else ""
            print(f"  {rank:02d}. {sig:<50} {stype:>5} | {change:.6f}{better_than_null}")

    # ----------------------------------------------------------------
    # 4. Cross-pair summary: average net movement per signature
    # ----------------------------------------------------------------
    print(f"\n{'='*80}")
    print("CROSS-PAIR SUMMARY (average movement toward P2 across all pairs)")
    print("="*80)

    all_sig_names = set()
    for pair_results in all_results.values():
        all_sig_names.update(pair_results.keys())

    cross_pair_avg = {}
    sig_type_map = {}
    for sig in all_sig_names:
        values = [
            all_results[pair][sig][0]
            for pair in all_results
            if sig in all_results[pair]
        ]
        # Collect sig_type from any pair that has it
        for pair in all_results:
            if sig in all_results[pair]:
                sig_type_map[sig] = all_results[pair][sig][1]
                break
        cross_pair_avg[sig] = (np.mean(values), np.std(values), len(values))

    sorted_cross = sorted(cross_pair_avg.items(), key=lambda x: x[1][0], reverse=True)
    null_cross_avg = cross_pair_avg.get('~~NULL (uniform random)', (0, 0, 0))[0]

    print(f"\n{'Rank':<5} {'Signature':<50} {'Type':>5} {'Avg Movement':>14} {'Std':>10} {'Pairs':>6} {'> Null?':>8}")
    print("-"*100)
    for rank, (sig, (avg, std, n)) in enumerate(sorted_cross, 1):
        better = "✓" if avg > null_cross_avg and sig != '~~NULL (uniform random)' else ""
        stype = sig_type_map.get(sig, '?')
        print(f"{rank:<5} {sig:<50} {stype:>5} {avg:>14.6f} {std:>10.6f} {n:>6}   {better}")

    print(f"\nNote: '✓' = outperforms null (uniform random) model on average.")
    print(f"Null model average movement: {null_cross_avg:.6f}")

    # ----------------------------------------------------------------
    # 5. Save all results to CSV
    # ----------------------------------------------------------------
    print(f"\n{'='*80}")
    print("SAVING RESULTS TO CSV")
    print("="*80)
    save_detailed_results_csv(all_results, output_dir)
    save_cross_pair_summary_csv(cross_pair_avg, sig_type_map, null_cross_avg, output_dir)
    save_wide_pivot_csv(all_results, output_dir)
    print(f"\nAll CSVs saved to: {output_dir}")
    print("  prism_detailed_results.csv   — one row per (pair, signature)")
    print("  prism_cross_pair_summary.csv — ranked leaderboard across all pairs")
    print("  prism_pivot_table.csv        — wide format, easy to open in Excel")

if __name__ == "__main__":
    main()