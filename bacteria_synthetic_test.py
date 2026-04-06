"""
prism_bacterial.py
==================
A PRISM-style protein mutation simulation using synthetic bacterial
mutational signatures instead of COSMIC cancer profiles.

Four synthetic SBS signatures are generated programmatically based on
known biophysical mechanisms:

  1. DMS  (Dimethyl Sulfate)
         — Alkylates bases; persistent AT bias.
           GC → AT transitions dominate (methylation blocks repair).
           P(G>A) and P(C>T) are heavily upweighted.

  2. Nitrous Acid
         — Deaminates cytosine to uracil (C→T) and adenine to
           hypoxanthine (A→G). Classic C→T and A→G transitions.
           Note: listed user bias is C→A, T→G which is modelled here.

  3. 5-Bromouracil (5-BrU)
         — Base analogue; mis-pairs as either A or G.
           Primary bias: C→A and T→G (same channel as nitrous acid
           but distinct mechanism; modelled with higher T→G weight).

  4. ROS  (Reactive Oxygen Species)
         — Oxidative burst. 8-oxoG causes G→T transversions.
           Strong GC→TA (G>T and C>A) with some G>C component.

Each synthetic signature is a dict { 'X[Y>Z]W': probability } in the
same trinucleotide SBS format used by the parent PRISM script, so the
existing apply_sbs_mutation() function works unchanged.

Usage
-----
  python prism_bacterial.py

The script is self-contained. Drop it in the same directory as your
PRISM project; it re-uses the helper functions from the original code
by importing them directly (see IMPORT block below).

If you prefer not to import from the original script, all required
helpers are duplicated at the bottom of this file under the comment
"# ---- Standalone fallback helpers ----".
"""

import os
import csv
import random
import numpy as np
import torch
from datetime import datetime
from scipy.spatial.distance import cosine

# -----------------------------------------------------------------------
# IMPORT from original PRISM script
# Adjust the path / module name to match your project layout.
# -----------------------------------------------------------------------
try:
    from prism import (
        load_esm_model,
        get_esm_embedding,
        apply_sbs_mutation,
        translate_nt_to_aa,
        mock_reverse_translate,
        read_fasta,
        apply_null_mutation,
    )
    print("Imported helpers from prism.py")
except ImportError:
    print("prism.py not found on path — using standalone fallbacks defined below.")
    # Fallback definitions are at the bottom of this file.
    # They will be registered into the namespace automatically.


# =======================================================================
# SYNTHETIC BACTERIAL MUTATIONAL SIGNATURES
# =======================================================================

# ------------------------------------------------------------------
# Internal helper: build a full 96-channel SBS profile from a small
# set of (ref_base, alt_base) → weight rules.
#
# The 96 SBS contexts are all combinations of:
#   left_base ∈ {A,C,G,T}  ×  ref ∈ {C,T}  ×  right_base ∈ {A,C,G,T}
# (COSMIC convention collapses 192 → 96 by always using the pyrimidine
# strand.  For biological simplicity here we use all 4 ref bases and
# all 4 left/right contexts.)
#
# Rule dict format:  { (ref_base, alt_base): relative_weight, ... }
# Unlisted pairs get a small background weight (epsilon) so the
# profile is never completely zero on any context.
# ------------------------------------------------------------------

BASES = list('ACGT')
EPSILON = 0.005   # background weight for unspecified substitutions


def _build_sbs_profile(rules: dict) -> dict:
    """
    Expand per-(ref,alt) weights into the full 192-channel trinucleotide
    SBS profile format  X[Y>Z]W  used by apply_sbs_mutation().

    rules: { (ref_base, alt_base): weight_float, ... }
    Returns: { 'A[C>A]A': prob, ... }   (normalised to sum=1)
    """
    raw = {}
    for left in BASES:
        for ref in BASES:
            for right in BASES:
                for alt in BASES:
                    if alt == ref:
                        continue
                    key = f"{left}[{ref}>{alt}]{right}"
                    weight = rules.get((ref, alt), EPSILON)
                    raw[key] = weight

    total = sum(raw.values())
    return {k: v / total for k, v in raw.items()}


# ------------------------------------------------------------------
# 1. DMS — Dimethyl Sulfate
#    Methylates N7-G and N3-A → stalls replication → GC→AT bias.
#    G>A  and  C>T  are the dominant transitions.
#    Low but non-zero rate for other changes (leaky repair).
# ------------------------------------------------------------------
DMS_RULES = {
    ('G', 'A'): 40.0,   # primary: methylated G mis-read as A
    ('C', 'T'): 35.0,   # complementary strand: G>A appears as C>T
    ('A', 'T'): 5.0,    # N3-A methylation → AT→TA transversion (minor)
    ('T', 'A'): 3.0,    # minor background
    ('G', 'T'): 2.0,    # minor transversion
    ('C', 'A'): 2.0,    # minor transversion
}

SIG_DMS = _build_sbs_profile(DMS_RULES)


# ------------------------------------------------------------------
# 2. Nitrous Acid
#    Deaminates C → U (read as T) and A → hypoxanthine (read as G).
#    User-specified bias: C→A and T→G.
# ------------------------------------------------------------------
NITROUS_RULES = {
    ('C', 'A'): 45.0,   # primary: deamination product mis-pair
    ('T', 'G'): 40.0,   # primary: A deamination on complementary strand
    ('G', 'T'): 4.0,    # minor
    ('A', 'C'): 2.0,    # minor
}

SIG_NITROUS = _build_sbs_profile(NITROUS_RULES)


# ------------------------------------------------------------------
# 3. 5-Bromouracil (5-BrU)
#    Base analogue that incorporates as T but can mis-pair with G.
#    Produces C→A and T→G transitions; T→G is slightly more prominent
#    than in nitrous acid because 5-BrU incorporates in place of T.
# ------------------------------------------------------------------
BROU_RULES = {
    ('C', 'A'): 30.0,   # 5-BrC mis-pair with A
    ('T', 'G'): 50.0,   # 5-BrU (incorporated as T) mis-pairs with G
    ('A', 'G'): 5.0,    # minor
    ('G', 'C'): 3.0,    # minor
}

SIG_BROU = _build_sbs_profile(BROU_RULES)


# ------------------------------------------------------------------
# 4. ROS — Reactive Oxygen Species (oxidative burst)
#    8-oxo-7,8-dihydroguanine (8-oxoG) pairs with A → G>T transversion.
#    Also produces C>A on the complementary strand.
#    Minor G>C component from strand-break repair.
# ------------------------------------------------------------------
ROS_RULES = {
    ('G', 'T'): 50.0,   # 8-oxoG : A mis-pair → G>T transversion
    ('C', 'A'): 45.0,   # complementary: C>A
    ('G', 'C'): 8.0,    # minor: 8-oxoG : C (correct) → strand-break artefact
    ('A', 'T'): 3.0,    # minor oxidative damage to A
}

SIG_ROS = _build_sbs_profile(ROS_RULES)


# ------------------------------------------------------------------
# Master registry of all synthetic signatures
# ------------------------------------------------------------------
BACTERIAL_SIGNATURES = {
    'DMS (GC→AT bias)':          (SIG_DMS,    'SBS'),
    'Nitrous Acid (C→A, T→G)':   (SIG_NITROUS,'SBS'),
    '5-BrU (C→A, T→G)':         (SIG_BROU,   'SBS'),
    'ROS burst (GC→TA)':         (SIG_ROS,    'SBS'),
}


# =======================================================================
# SIMULATION RUNNER  (mirrors run_simulation() in prism.py)
# =======================================================================

def run_simulation_bacterial(p1_nt_orig, p2_embedding, p1_orig_embedding,
                              profile, model, alphabet, batch_converter,
                              max_steps, use_null=False):
    """
    Iteratively mutates P1's nucleotide sequence using the given SBS profile
    (or uniform random if use_null=True) and tracks cosine distance to P2.

    Returns: list of cosine distances at each step (index 0 = initial).
    """
    initial_distance = cosine(p1_orig_embedding, p2_embedding)
    current_nt = p1_nt_orig
    trajectory = [initial_distance]

    for step in range(1, max_steps + 1):
        if use_null:
            current_nt, mutated = apply_null_mutation(current_nt)
        else:
            current_nt, mutated = apply_sbs_mutation(current_nt, profile)

        if not mutated:
            print("  -> Halted early: no valid mutation contexts remain.")
            break

        current_aa = translate_nt_to_aa(current_nt)
        if not current_aa:
            break

        current_emb = get_esm_embedding(current_aa, model, alphabet, batch_converter)
        trajectory.append(cosine(current_emb, p2_embedding))

    return trajectory


# =======================================================================
# CSV SAVERS
# =======================================================================

def save_bacterial_detailed_csv(all_results, output_dir):
    path = os.path.join(output_dir, 'bacterial_detailed_results.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['protein_pair', 'signature', 'sig_type', 'net_movement_toward_p2'])
        for pair, pair_res in all_results.items():
            for sig, (net, stype) in pair_res.items():
                w.writerow([pair, sig, stype, f"{net:.6f}"])
    print(f"  Saved detailed results  -> {path}")
    return path


def save_bacterial_summary_csv(cross_avg, sig_type_map, null_avg, output_dir):
    path = os.path.join(output_dir, 'bacterial_cross_pair_summary.csv')
    sorted_cross = sorted(cross_avg.items(), key=lambda x: x[1][0], reverse=True)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank','signature','sig_type','avg_movement_toward_p2',
                    'std','n_pairs','beats_null_model'])
        for rank, (sig, (avg, std, n)) in enumerate(sorted_cross, 1):
            beats = 'YES' if avg > null_avg and sig != '~~NULL (uniform random)' else 'NO'
            w.writerow([rank, sig, sig_type_map.get(sig,'?'),
                        f"{avg:.6f}", f"{std:.6f}", n, beats])
    print(f"  Saved cross-pair summary -> {path}")
    return path


def save_bacterial_pivot_csv(all_results, output_dir):
    path = os.path.join(output_dir, 'bacterial_pivot_table.csv')
    pair_labels = list(all_results.keys())
    all_sigs = sorted({sig for pr in all_results.values() for sig in pr})
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['signature'] + pair_labels + ['avg_across_pairs'])
        for sig in all_sigs:
            vals = [all_results[p][sig][0] if sig in all_results[p] else '' for p in pair_labels]
            nums = [v for v in vals if v != '']
            avg  = np.mean(nums) if nums else ''
            w.writerow([sig] + [f"{v:.6f}" if v!='' else '' for v in vals]
                       + [f"{avg:.6f}" if avg!='' else ''])
    print(f"  Saved pivot table        -> {path}")
    return path


# =======================================================================
# MAIN
# =======================================================================

def main():
    base_path = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM'
    MAX_STEPS    = 50
    NULL_REPS    = 5   # replicates for null model averaging

    timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(base_path, f'bacterial_results_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Results will be saved to: {output_dir}")

    # ----------------------------------------------------------------
    # Protein pairs — same format as prism.py
    # ----------------------------------------------------------------
    PROTEIN_PAIRS = {
        'PAIR1-PDK': ('A0A010RJQ9.fa', 'A0A010S4M5.fa', 'PAIR1-PDK'),
        'PAIR2-Ras': ('A0A010R5C7.fa', 'A0A010SAT4.fa', 'PAIR2-Ras'),
        'PAIR3-p53': ('A0A087VQ67.fa', 'A0A087WT22.fa', 'PAIR3-p53'),
        'PAIR4-ABC': ('A0A010Z411.fa', 'A0A011NXC8.fa', 'PAIR4-ABC'),
        'PAIR5-HR':  ('A0A009T0V7.fa', 'A0A010QZB6.fa', 'PAIR5-HR'),
    }

    # ----------------------------------------------------------------
    # Print synthetic signature summaries
    # ----------------------------------------------------------------
    print("\n" + "="*70)
    print("SYNTHETIC BACTERIAL MUTATIONAL SIGNATURES")
    print("="*70)
    for name, (profile, stype) in BACTERIAL_SIGNATURES.items():
        # Show top-5 highest-probability contexts as a sanity check
        top5 = sorted(profile.items(), key=lambda x: x[1], reverse=True)[:5]
        print(f"\n  {name}  [{stype}]")
        for ctx, p in top5:
            print(f"    {ctx}  →  {p:.4f}")

    # ----------------------------------------------------------------
    # Load ESM-2 once
    # ----------------------------------------------------------------
    model, alphabet, batch_converter = load_esm_model()

    all_results = {}

    # ----------------------------------------------------------------
    # Loop over protein pairs
    # ----------------------------------------------------------------
    for pair_label, (p1_file, p2_file, subfolder) in PROTEIN_PAIRS.items():
        print(f"\n{'='*70}")
        print(f"PROTEIN PAIR: {pair_label}")
        print(f"{'='*70}")

        p1_path = os.path.join(base_path, 'protein_pairs', subfolder, p1_file)
        p2_path = os.path.join(base_path, 'protein_pairs', subfolder, p2_file)

        if not os.path.exists(p1_path) or not os.path.exists(p2_path):
            print(f"  WARNING: Files not found for {pair_label}, skipping.")
            continue

        p1_aa_orig = read_fasta(p1_path)
        p2_aa      = read_fasta(p2_path)

        p1_nt_orig       = mock_reverse_translate(p1_aa_orig)
        p1_aa_sync       = translate_nt_to_aa(p1_nt_orig)

        p2_emb           = get_esm_embedding(p2_aa,      model, alphabet, batch_converter)
        p1_orig_emb      = get_esm_embedding(p1_aa_sync, model, alphabet, batch_converter)
        initial_distance = cosine(p1_orig_emb, p2_emb)

        print(f"  Initial cosine distance (P1 → P2): {initial_distance:.6f}")

        pair_results = {}

        # ---- A. Bacterial signature runs ----
        for sig_name, (profile, sig_type) in BACTERIAL_SIGNATURES.items():
            print(f"\n  --- Signature: {sig_name} ---")
            traj = run_simulation_bacterial(
                p1_nt_orig, p2_emb, p1_orig_emb,
                profile, model, alphabet, batch_converter,
                MAX_STEPS, use_null=False
            )
            net = initial_distance - traj[-1]
            pair_results[sig_name] = (net, sig_type)
            print(f"  Final distance: {traj[-1]:.6f} | Movement toward P2: {net:.6f}")

        # ---- B. Null model ----
        print(f"\n  --- NULL MODEL (uniform, {NULL_REPS} replicates) ---")
        null_vals = []
        for _ in range(NULL_REPS):
            traj = run_simulation_bacterial(
                p1_nt_orig, p2_emb, p1_orig_emb,
                profile=None, model=model, alphabet=alphabet,
                batch_converter=batch_converter,
                max_steps=MAX_STEPS, use_null=True
            )
            null_vals.append(initial_distance - traj[-1])

        null_avg = np.mean(null_vals)
        null_std = np.std(null_vals)
        pair_results['~~NULL (uniform random)'] = (null_avg, 'NULL')
        print(f"  Null avg movement: {null_avg:.6f} ± {null_std:.6f}")

        all_results[pair_label] = pair_results

        # Per-pair leaderboard
        print(f"\n  LEADERBOARD — {pair_label}")
        print(f"  {'Signature':<45} {'Type':>5} | Net movement toward P2")
        print(f"  {'-'*65}")
        for rank, (sig, (chg, stype)) in enumerate(
                sorted(pair_results.items(), key=lambda x: x[1][0], reverse=True), 1):
            marker = " ✓" if chg > null_avg and sig != '~~NULL (uniform random)' else ""
            print(f"  {rank:02d}. {sig:<45} {stype:>5} | {chg:.6f}{marker}")

    # ----------------------------------------------------------------
    # Cross-pair summary
    # ----------------------------------------------------------------
    print(f"\n{'='*80}")
    print("CROSS-PAIR SUMMARY")
    print("="*80)

    all_sig_names = {sig for pr in all_results.values() for sig in pr}
    cross_avg    = {}
    sig_type_map = {}

    for sig in all_sig_names:
        vals = [all_results[p][sig][0] for p in all_results if sig in all_results[p]]
        cross_avg[sig] = (np.mean(vals), np.std(vals), len(vals))
        for p in all_results:
            if sig in all_results[p]:
                sig_type_map[sig] = all_results[p][sig][1]
                break

    null_cross_avg = cross_avg.get('~~NULL (uniform random)', (0, 0, 0))[0]
    sorted_cross   = sorted(cross_avg.items(), key=lambda x: x[1][0], reverse=True)

    print(f"\n{'Rank':<5} {'Signature':<45} {'Type':>5} "
          f"{'Avg Movement':>14} {'Std':>10} {'Pairs':>6} {'> Null?':>8}")
    print("-"*95)
    for rank, (sig, (avg, std, n)) in enumerate(sorted_cross, 1):
        better = "✓" if avg > null_cross_avg and sig != '~~NULL (uniform random)' else ""
        print(f"{rank:<5} {sig:<45} {sig_type_map.get(sig,'?'):>5} "
              f"{avg:>14.6f} {std:>10.6f} {n:>6}   {better}")

    print(f"\nNull model average movement: {null_cross_avg:.6f}")

    # ----------------------------------------------------------------
    # Save CSVs
    # ----------------------------------------------------------------
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print("="*80)
    save_bacterial_detailed_csv(all_results, output_dir)
    save_bacterial_summary_csv(cross_avg, sig_type_map, null_cross_avg, output_dir)
    save_bacterial_pivot_csv(all_results, output_dir)
    print(f"\nAll CSVs saved to: {output_dir}")


# =======================================================================
# STANDALONE FALLBACK HELPERS
# (used automatically if prism.py cannot be imported)
# =======================================================================

try:
    # Already imported above — skip redefinition
    apply_null_mutation
except NameError:
    import esm
    from Bio.Seq import Seq

    CODON_MAP = {
        'A':'GCC','C':'TGC','D':'GAC','E':'GAG','F':'TTC',
        'G':'GGC','H':'CAC','I':'ATC','K':'AAG','L':'CTG',
        'M':'ATG','N':'AAC','P':'CCC','Q':'CAG','R':'CGC',
        'S':'AGC','T':'ACC','V':'GTG','W':'TGG','Y':'TAC','*':'TAA'
    }

    def mock_reverse_translate(aa_seq):
        return "".join([CODON_MAP.get(aa,'NNN')
                        for aa in aa_seq.upper().replace('\n','') if aa in CODON_MAP])

    def load_esm_model():
        print("Loading ESM-2 model...")
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
        reps = results["representations"][33]
        return reps[0, 1:len(aa_seq)+1].mean(0).numpy()

    def apply_sbs_mutation(nt_seq, profile):
        seq_list = list(nt_seq.upper())
        possible = []
        for idx in range(len(seq_list) - 2):
            ctx = "".join(seq_list[idx:idx+3])
            for mut_sig, prob in profile.items():
                if len(mut_sig) < 7 or prob <= 0:
                    continue
                orig_tn = mut_sig[0] + mut_sig[2] + mut_sig[6]
                if ctx == orig_tn:
                    possible.append({'index':idx+1,'new_base':mut_sig[4],'weight':prob})
        if not possible:
            return nt_seq, False
        weights = [m['weight'] for m in possible]
        total = sum(weights)
        probs = [w/total for w in weights]
        chosen = np.random.choice(possible, p=probs)
        seq_list[chosen['index']] = chosen['new_base']
        return "".join(seq_list), True

    def apply_null_mutation(nt_seq):
        seq_list = list(nt_seq.upper())
        idx = random.randint(0, len(seq_list)-1)
        original = seq_list[idx]
        bases = [b for b in 'ACGT' if b != original]
        seq_list[idx] = random.choice(bases)
        return "".join(seq_list), True

    def translate_nt_to_aa(nt_seq):
        trimmed = nt_seq[:len(nt_seq)-(len(nt_seq)%3)]
        return str(Seq(trimmed).translate(to_stop=True))

    def read_fasta(filepath):
        with open(filepath,'r') as f:
            lines = f.readlines()
        return "".join([l.strip() for l in lines if not l.startswith('>')])


if __name__ == "__main__":
    main()