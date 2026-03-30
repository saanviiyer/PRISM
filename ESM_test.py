import os
import random
from Bio.Seq import Seq
import torch
import esm
from pathlib import Path

# ---- ESM setup (load once, reuse) ----
def load_esm_model():
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    return model, alphabet, batch_converter

def esm_score_mutations(wt_aa_seq, mut_aa_seq, model, alphabet, batch_converter):
    """
    Compute masked marginal log-likelihood ratio for each mutated position.
    Returns a dict of {position: {wt, mut, llr}} and a mean score.
    Negative LLR = mutation is evolutionarily unusual (likely deleterious).
    """
    wt_seq = wt_aa_seq.upper()
    mut_seq = mut_aa_seq.upper()

    changed_positions = [i for i, (a, b) in enumerate(zip(wt_seq, mut_seq)) if a != b]
    if not changed_positions:
        return {}, 0.0

    llr_scores = {}

    for pos in changed_positions:
        wt_aa = wt_seq[pos]
        mut_aa = mut_seq[pos]

        masked = list(wt_seq)
        masked[pos] = alphabet.mask_idx

        masked_seq = wt_seq[:pos] + "<mask>" + wt_seq[pos+1:]

        data = [("wt", masked_seq)]
        _, _, tokens = batch_converter(data)

        with torch.no_grad():
            results = model(tokens, repr_layers=[], return_contacts=False)

        logits = results["logits"][0, pos+1]
        log_probs = torch.log_softmax(logits, dim=-1)

        wt_idx = alphabet.get_idx(wt_aa)
        mut_idx = alphabet.get_idx(mut_aa)

        llr = (log_probs[mut_idx] - log_probs[wt_idx]).item()
        # Store the specific change alongside the score
        llr_scores[pos] = {"wt": wt_aa, "mut": mut_aa, "llr": llr}

    mean_llr = sum(m["llr"] for m in llr_scores.values()) / len(llr_scores)
    return llr_scores, mean_llr


# ---- Translation helper ----
def translate_nt_to_aa(nt_seq):
    trimmed = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)]
    aa_seq = str(Seq(trimmed).translate(to_stop=True))
    return aa_seq


# ---- Pipeline ----
def get_mutation_rate_with_functional_scoring(
    protein_nt_seq,
    file_path,
    esm_model=None,
    alphabet=None,
    batch_converter=None
):
    mutation_signatures = {}
    with open(file_path, 'r') as f:
        next(f, None)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                mutation_signatures[parts[0]] = float(parts[1])

    seq_list = list(protein_nt_seq.upper())
    mutations_count = 0

    for idx in range(len(seq_list) - 2):
        context = "".join(seq_list[idx:idx+3])
        for mut_sig, prob in mutation_signatures.items():
            orig_tn = mut_sig[0] + mut_sig[2] + mut_sig[6]
            if context == orig_tn:
                if random.random() < prob:
                    new_base = mut_sig[4]
                    seq_list[idx+1] = new_base
                    mutations_count += 1
                    break

    mutated_nt_seq = "".join(seq_list)
    mutation_rate = mutations_count / len(protein_nt_seq) if protein_nt_seq else 0

    esm_mean_llr = None
    detailed_llrs = {}
    if esm_model is not None:
        wt_aa = translate_nt_to_aa(protein_nt_seq)
        mut_aa = translate_nt_to_aa(mutated_nt_seq)
        detailed_llrs, esm_mean_llr = esm_score_mutations(
            wt_aa, mut_aa, esm_model, alphabet, batch_converter
        )

    return mutation_rate, esm_mean_llr, detailed_llrs


# ---- Main analysis loop ----
profile_directory = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/profiles'
protein_directory  = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/proteins'
USE_ESM = True

categories = ['Tobacco', 'UV', 'Mismatch', 'Chemotherapy']
mut_rates = {}

# Global tracker for the most damaging mutations across the entire run
global_worst_mutations = []

esm_model, alphabet, batch_converter = (None, None, None)
if USE_ESM:
    print("Loading ESM-2 model...")
    esm_model, alphabet, batch_converter = load_esm_model()

print("Starting mutation simulation analysis...\n" + "-"*40)

for category in categories:
    cat_profile_dir = os.path.join(profile_directory, category)
    cat_protein_dir = os.path.join(protein_directory, category)

    if not os.path.exists(cat_profile_dir) or not os.path.exists(cat_protein_dir):
        continue

    profiles = [f for f in os.listdir(cat_profile_dir) if not f.startswith('.')]
    proteins = [f for f in os.listdir(cat_protein_dir) if not f.startswith('.')]

    for profile in profiles:
        profile_path = os.path.join(cat_profile_dir, profile)

        for protein_file in proteins:
            protein_path = os.path.join(cat_protein_dir, protein_file)

            with open(protein_path, 'r') as f:
                lines = f.readlines()
                protein_seq = "".join([l.strip() for l in lines if not l.startswith('>')])

            # Unpack the new detailed_llrs dictionary
            rate, llr, detailed_llrs = get_mutation_rate_with_functional_scoring(
                protein_seq,
                profile_path,
                esm_model=esm_model,
                alphabet=alphabet,
                batch_converter=batch_converter
            )

            result_key = f"{category.upper()} | Profile: {profile} | Protein: {protein_file}"
            llr_str = f"{llr:.4f}" if llr is not None else "N/A"
            
            print(f"\n{result_key}")
            print(f"  Overall Rate: {rate:.6f} | Mean ESM LLR: {llr_str}")

            # Identify and print the most damaging mutations for THIS protein
            if detailed_llrs:
                # Sort mutations by LLR score (ascending: most negative/damaging first)
                sorted_muts = sorted(detailed_llrs.items(), key=lambda x: x[1]["llr"])
                
                print("  Top 3 Most Damaging Mutations:")
                for pos, data in sorted_muts[:3]:
                    print(f"    -> Pos {pos}: {data['wt']} changed to {data['mut']} (LLR: {data['llr']:.4f})")

                # Add to the global tracking list
                for pos, data in detailed_llrs.items():
                    global_worst_mutations.append({
                        "category": category,
                        "protein": protein_file,
                        "pos": pos,
                        "wt": data["wt"],
                        "mut": data["mut"],
                        "llr": data["llr"]
                    })

# ---- Final Global Output ----
print("\n" + "="*50)
print(" TOP 10 MOST DAMAGING MUTATIONS OVERALL")
print("="*50)

# Sort the global list by LLR (most negative first)
global_worst_mutations.sort(key=lambda x: x["llr"])

for i, mut in enumerate(global_worst_mutations[:10], 1):
    print(f"{i}. Protein: {mut['protein']} ({mut['category']})")
    print(f"   Mutation: Pos {mut['pos']} [{mut['wt']} -> {mut['mut']}] | LLR: {mut['llr']:.4f}\n")