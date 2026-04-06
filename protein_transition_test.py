import os
import random
import torch
import esm
import numpy as np
from Bio.Seq import Seq
from scipy.spatial.distance import cosine

# ---- Mock Reverse Translation (Scaffolding for Testing) ----
# Maps Amino Acids to their most common human codons
CODON_MAP = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGC',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
    '*': 'TAA'
}

def mock_reverse_translate(aa_seq):
    """Converts Amino Acid sequence to a plausible DNA sequence."""
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
    aa_seq = aa_seq[:1022] # Truncate to ESM max length
    data = [("seq", aa_seq.upper())]
    _, _, tokens = batch_converter(data)
    with torch.no_grad():
        results = model(tokens, repr_layers=[33], return_contacts=False)
    token_representations = results["representations"][33]
    return token_representations[0, 1 : len(aa_seq) + 1].mean(0).numpy()

# ---- Dynamic Profile Discovery ----
def discover_profiles(base_dir):
    """Automatically finds all .txt signature files in the directory tree."""
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
def load_mutation_profile(profile_path):
    sigs = {}
    with open(profile_path, 'r') as f:
        next(f, None)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sigs[parts[0]] = float(parts[1])
    return sigs

def apply_single_weighted_mutation(nt_seq, profile):
    seq_list = list(nt_seq.upper())
    possible_mutations = []
    
    for idx in range(len(seq_list) - 2):
        context = "".join(seq_list[idx:idx+3])
        for mut_sig, prob in profile.items():
            orig_tn = mut_sig[0] + mut_sig[2] + mut_sig[6]
            if context == orig_tn and prob > 0:
                possible_mutations.append({
                    'index': idx + 1,
                    'new_base': mut_sig[4],
                    'weight': prob
                })
                
    if not possible_mutations:
        return nt_seq, False

    weights = [m['weight'] for m in possible_mutations]
    probs = [w / sum(weights) for w in weights]
    chosen_mut = np.random.choice(possible_mutations, p=probs)
    
    seq_list[chosen_mut['index']] = chosen_mut['new_base']
    return "".join(seq_list), True

def translate_nt_to_aa(nt_seq):
    trimmed = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)]
    return str(Seq(trimmed).translate(to_stop=True))

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        return "".join([l.strip() for l in lines if not l.startswith('>')])

# ---- Main Simulation Loop ----
def main():
    # Absolute paths based on your environment
    base_path = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM'
    
    p1_aa_path = os.path.join(base_path, 'protein_pairs/PAIR1/A0A010RJQ9.fa')
    p2_aa_path = os.path.join(base_path, 'protein_pairs/PAIR1/A0A010S4M5.fa')
    profile_base_directory = os.path.join(base_path, 'profiles')

    MAX_STEPS = 50

    # 1. Dynamically load all profiles
    PROFILES_TO_TEST = discover_profiles(profile_base_directory)
    if not PROFILES_TO_TEST:
        print("No profiles found to test. Check your 'profiles' folder path.")
        return
    print(f"Found {len(PROFILES_TO_TEST)} mutational profiles to test.")

    # 2. Load Models
    model, alphabet, batch_converter = load_esm_model()
    
    # 3. Load Sequences and Apply Mock Translation to P1
    p1_aa_orig = read_fasta(p1_aa_path)
    p2_aa = read_fasta(p2_aa_path)
    
    print("\nApplying mock reverse-translation to P1...")
    p1_nt_orig = mock_reverse_translate(p1_aa_orig)
    
    # Recalculate P1 AA from the new DNA to ensure perfect syncing
    p1_aa_sync = translate_nt_to_aa(p1_nt_orig)
    
    p2_embedding = get_esm_embedding(p2_aa, model, alphabet, batch_converter)
    p1_orig_embedding = get_esm_embedding(p1_aa_sync, model, alphabet, batch_converter)
    initial_distance = cosine(p1_orig_embedding, p2_embedding)

    print("\n" + "="*60)
    print(f"STARTING COMPARATIVE SIMULATION")
    print(f"Initial distance (P1 -> P2): {initial_distance:.6f}")
    print("="*60)

    results_summary = {}

    # 4. Test every discovered signature
    for signature_name, profile_path in sorted(PROFILES_TO_TEST.items()):
        print(f"\n--- Running: {signature_name} ---")
        profile = load_mutation_profile(profile_path)
        
        current_nt = p1_nt_orig 
        trajectory = [initial_distance]
        
        for step in range(1, MAX_STEPS + 1):
            current_nt, mutated = apply_single_weighted_mutation(current_nt, profile)
            if not mutated:
                print("  -> Simulation halted early: No valid context mutations left.")
                break
                
            current_aa = translate_nt_to_aa(current_nt)
            current_embedding = get_esm_embedding(current_aa, model, alphabet, batch_converter)
            
            dist = cosine(current_embedding, p2_embedding)
            trajectory.append(dist)
        
        net_change = initial_distance - trajectory[-1]
        results_summary[signature_name] = net_change
        
        print(f"Final Distance after {len(trajectory)-1} mutations: {trajectory[-1]:.6f}")
        print(f"Total movement toward P2: {net_change:.6f}")

    # 5. Print the final leaderboard
    print("\n" + "="*80)
    print("EFFICIENCY LEADERBOARD (Which signature moved P1 closest to P2?)")
    print("="*80)
    sorted_results = sorted(results_summary.items(), key=lambda x: x[1], reverse=True)
    for rank, (sig, change) in enumerate(sorted_results, 1):
        print(f"{rank:02d}. {sig:<40} | Moved {change:.6f} closer")

if __name__ == "__main__":
    main()