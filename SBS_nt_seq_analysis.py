# 3/23 AITHYRA PROJECT PROOF OF CONCEPT
import webbrowser
import random


# Open the URL in the default web browser
# sample protein: YTHDF2 (a gatekeeper protein preventing UV-induced inflammation)

# YTHDF2:str = """GCCGCGCGCTGTGTCTCCGCTGCGTCCGCCGAGGCCCCCGAGTGTCAGGGACAAAAGCCT
# CCGCCTGCTCCCGCAGCCGGGGCTCATCTGCCGCCGCCGCCGCGCTGAGGAGAGTTCGCC
# GCCGTCGCCGCCCGTGAGGATCTGAGAGCCATGTCGGCCAGCAGCCTCTTGGAGCAGAGA
# CCAAAAGGTCAAGGAAACAAAGTACAAAATGGATCTGTACATCAAAAGGATGGATTAAAC
# GATGATGATTTTGAACCTTACTTGAGTCCACAGGCAAGGCCCAATAATGCATATACTGCC
# ATGTCAGATTCCTACTTACCCAGTTACTACAGTCCCTCCATTGGCTTCTCCTATTCTTTG
# GGTGAAGCTGCTTGGTCTACGGGGGGTGACACAGCCATGCCCTACTTAACTTCTTATGGA
# CAGCTGAGCAACGGAGAGCCCCACTTCCTACCAGATGCAATGTTTGGGCAACCAGGAGCC
# CTAGGTAGCACTCCATTTCTTGGTCAGCATGGTTTTAATTTCTTTCCCAGTGGGATTGAC
# TTCTCAGCATGGGGAAATAACAGTTCTCAGGGACAGTCTACTCAGAGCTCTGGATATAGT
# AGCAATTATGCTTATGCACCTAGCTCCTTAGGTGGAGCCATGATTGATGGACAGTCAGCT
# TTTGCCAATGAGACCCTCAATAAGGCTCCTGGCATGAATACTATAGACCAAGGGATGGCA
# GCACTGAAGTTGGGTAGCACAGAAGTTGCAAGCAATGTTCCAAAAGTTGTAGGTTCTGCT
# GTTGGTAGCGGGTCCATTACTAGTAACATCGTGGCTTCCAATAGTTTGCCTCCAGCCACC
# ATTGCTCCTCCAAAACCAGCATCTTGGGCTGATATTGCTAGCAAGCCTGCAAAACAGCAA
# CCTAAACTGAAGACCAAGAATGGCATTGCAGGGTCAAGTCTTCCGCCACCCCCGATAAAG
# CATAACATGGATATTGGAACTTGGGATAACAAGGGTCCCGTTGCAAAAGCCCCCTCACAG
# GCTTTGGTTCAGAATATAGGTCAGCCAACCCAGGGGTCTCCTCAGCCTGTAGGTCAGCAG
# GCTAACAATAGCCCACCAGTGGCTCAGGCATCAGTAGGGCAACAGACACAGCCATTGCCT
# CCACCTCCACCACAGCCTGCCCAGCTTTCAGTCCAGCAACAGGCAGCTCAGCCAACCCGC
# TGGGTAGCACCTCGGAACCGTGGCAGTGGGTTCGGTCATAATGGGGTGGATGGTAATGGA
# GTAGGACAGTCTCAGGCTGGTTCTGGATCTACTCCTTCAGAACCCCACCCAGTGTTGGAG
# AAGCTTCGGTCCATTAATAACTATAACCCCAAAGATTTTGACTGGAATCTGAAACATGGC
# CGGGTTTTCATCATTAAGAGCTACTCTGAGGACGATATTCACCGTTCCATTAAGTATAAT
# ATTTGGTGCAGCACAGAGCATGGTAACAAGAGACTGGATGCTGCTTATCGTTCCATGAAC
# GGGAAAGGCCCCGTTTACTTACTTTTCAGTGTCAACGGCAGTGGACACTTCTGTGGCGTG
# GCAGAAATGAAATCTGCTGTGGACTACAACACATGTGCAGGTGTGTGGTCCCAGGACAAA
# TGGAAGGGTCGTTTTGATGTCAGGTGGATTTTTGTGAAGGACGTTCCCAATAGCCAACTG
# CGACACATTCGCCTAGAGAACAACGAGAATAAACCAGTGACCAACTCTAGGGACACTCAG
# GAAGTGCCTCTGGAAAAGGCTAAGCAGGTGTTGAAAATTATAGCCAGCTACAAGCACACC
# ACTTCCATTTTTGATGACTTCTCACACTATGAGAAACGCCAAAGAGGAAGAAGAAAGTGT
# TAAAAAGGAACGTCAAGGTCGTGGGAAATAAAAGGCAGTTCTACACAGACTGCAGCAACG
# GTTGCATCTGCATATCCTAAGAGGAAAAAATGACCTTCAAGAGAATTAGGACTTTTTTCT
# TAATTTCACTGACTTCAGAGACGATTGCAGACTTGCAGTTTAAGTATTGGAATTTCACAA
# AAGACATAGGACTTAACTGGAAAATGAAAA"""

# # apply SBS 7a - ultraviolet light exposure
# mutation_signatures = {}

# file_path = 'SBS7a_PROFILE.txt'
# with open(file_path, 'r') as f:
#     for line in f:
#         if line.strip().split()[0] != 'SBS7a':
#             mutation_signatures[line.strip().split()[0]] = line.strip().split()[1]
# mutation_signatures.pop(next(iter(mutation_signatures)))

# # {'A[C>A]A': '6.70435065427006e-05', 'A[C>A]C': '0.000179116233897663', 'A[C>A]G': '7.12462338185117e-05', 'A[C>A]T': '0.000248161039143131', 'A[C>G]A': '6.49421429047951e-05', 'A[C>G]C': '2.16140259898856e-05', 'A[C>G]G': '5.08329870502865e-05', 'A[C>G]T': '9.56620779922714e-05', 'A[C>T]A': '0.000122079220868798', 'A[C>T]C': '0.0177114935194896', 'A[C>T]G': '2.23144805358541e-16', 'A[C>T]T': '0.00742481818726625', 'A[T>A]A': '0.00194125974168417', 'A[T>A]C': '0.000670435065427006', 'A[T>A]G': '0.000587381169262168', 'A[T>A]T': '0.00371240909363312', 'A[T>C]A': '0.00049031818217796', 'A[T>C]C': '0.000227147402764075', 'A[T>C]G': '0.000187121428708732', 'A[T>C]T': '0.000615399351100909', 'A[T>G]A': '0.000473307143204439', 'A[T>G]C': '0.000247160389791747', 'A[T>G]G': '0.000276179220981871', 'A[T>G]T': '0.00030920064957753', 'C[C>A]A': '0.000455295454879534', 'C[C>A]C': '0.000240155844332062', 'C[C>A]G': '2.23144805358541e-16', 'C[C>A]T': '0.000372241558714696', 'C[C>G]A': '2.23144805358541e-16', 'C[C>G]C': '0.000131085065031251', 'C[C>G]G': '7.75503247322284e-05', 'C[C>G]T': '9.56620779922714e-05', 'C[C>T]A': '0.0500324675691796', 'C[C>T]C': '0.0553359091315126', 'C[C>T]G': '0.0067143571477839', 'C[C>T]T': '0.0566367532883113', 'C[T>A]A': '0.000177114935194896', 'C[T>A]C': '0.000261169480711117', 'C[T>A]G': '0.000112072727354962', 'C[T>A]T': '0.000371240909363312', 'C[T>C]A': '0.000451292857474', 'C[T>C]C': '0.00025916818200835', 'C[T>C]G': '0.000104067532543894', 'C[T>C]T': '0.00182118181951814', 'C[T>G]A': '8.42546753864984e-05', 'C[T>G]C': '0.000172111688437978', 'C[T>G]G': '2.23144805358541e-16', 'C[T>G]T': '2.65172078116652e-05', 'G[C>A]A': '4.06263636661738e-05', 'G[C>A]C': '0.000138089610490936', 'G[C>A]G': '4.25275974338026e-05', 'G[C>A]T': '2.23144805358541e-16', 'G[C>G]A': '1.96127272871184e-05', 'G[C>G]C': '3.32215584659352e-05', 'G[C>G]G': '4.51292857474e-05', 'G[C>G]T': '4.20272727581108e-05', 'G[C>T]A': '8.49551299324669e-05', 'G[C>T]C': '0.011707597411188', 'G[C>T]G': '5.52358441963743e-06', 'G[C>T]T': '0.0161104545572758', 'G[T>A]A': '0.00040126038990482', 'G[T>A]C': '0.000224145454709925', 'G[T>A]G': '0.000290188311901242', 'G[T>A]T': '0.000306198701523379', 'G[T>C]A': '0.000101065584489743', 'G[T>C]C': '7.41481169375241e-05', 'G[T>C]G': '7.25470779753104e-05', 'G[T>C]T': '0.000811526623972093', 'G[T>G]A': '9.06588312353534e-05', 'G[T>G]C': '0.000121078571517415', 'G[T>G]G': '0.000325211039199667', 'G[T>G]T': '0.000835542208405299', 'T[C>A]A': '0.000390253247039601', 'T[C>A]C': '0.000636412987479964', 'T[C>A]G': '3.44223376875956e-05', 'T[C>A]T': '0.000613398052398142', 'T[C>G]A': '2.23144805358541e-16', 'T[C>G]C': '5.09330519854248e-05', 'T[C>G]G': '6.35412338128581e-06', 'T[C>G]T': '0.000281182467738789', 'T[C>T]A': '0.238154545629295', 'T[C>T]C': '0.331214935307969', 'T[C>T]G': '0.0741481169375241', 'T[C>T]T': '0.107069480598044', 'T[T>A]A': '0.000356231169092559', 'T[T>A]C': '0.000389252597688217', 'T[T>A]G': '2.23144805358541e-16', 'T[T>A]T': '0.000798518182404106', 'T[T>C]A': '0.00215139610547472', 'T[T>C]C': '0.000596387013424621', 'T[T>C]G': '0.000319207143091366', 'T[T>C]T': '0.000812527273323476', 'T[T>G]A': '0.0001280831169771', 'T[T>G]C': '0.000116075324760497', 'T[T>G]G': '2.23144805358541e-16', 'T[T>G]T': '8.29538312296997e-05'}

# YTHDF2_mut = YTHDF2

# for i in mutation_signatures.keys():
#     orig_tn = i[0] + i[2] + i[6]
#     new_tn = i[0] + i[4] + i[6]

#     prob = float(mutation_signatures[i])
    
#     if random.random() < (prob * 100):
#         YTHDF2 = YTHDF2_mut.replace(orig_tn, new_tn)

# print(YTHDF2)

# count = sum(1 for a, b in zip(YTHDF2, YTHDF2_mut) if a != b)

# mutation_rate = count / len(YTHDF2)
# print(mutation_rate)

# def get_mutation_rate(protein_nt_seq, file_path):
#     mutation_signatures = {}
    
#     with open(file_path, 'r') as f:
    
#         for line in f:
#             if line.strip().split()[0] != 'SBS7a':
#                 mutation_signatures[line.strip().split()[0]] = line.strip().split()[1]
#     mutation_signatures.pop(next(iter(mutation_signatures)))

#     # {'A[C>A]A': '6.70435065427006e-05', 'A[C>A]C': '0.000179116233897663', 'A[C>A]G': '7.12462338185117e-05', 'A[C>A]T': '0.000248161039143131', 'A[C>G]A': '6.49421429047951e-05', 'A[C>G]C': '2.16140259898856e-05', 'A[C>G]G': '5.08329870502865e-05', 'A[C>G]T': '9.56620779922714e-05', 'A[C>T]A': '0.000122079220868798', 'A[C>T]C': '0.0177114935194896', 'A[C>T]G': '2.23144805358541e-16', 'A[C>T]T': '0.00742481818726625', 'A[T>A]A': '0.00194125974168417', 'A[T>A]C': '0.000670435065427006', 'A[T>A]G': '0.000587381169262168', 'A[T>A]T': '0.00371240909363312', 'A[T>C]A': '0.00049031818217796', 'A[T>C]C': '0.000227147402764075', 'A[T>C]G': '0.000187121428708732', 'A[T>C]T': '0.000615399351100909', 'A[T>G]A': '0.000473307143204439', 'A[T>G]C': '0.000247160389791747', 'A[T>G]G': '0.000276179220981871', 'A[T>G]T': '0.00030920064957753', 'C[C>A]A': '0.000455295454879534', 'C[C>A]C': '0.000240155844332062', 'C[C>A]G': '2.23144805358541e-16', 'C[C>A]T': '0.000372241558714696', 'C[C>G]A': '2.23144805358541e-16', 'C[C>G]C': '0.000131085065031251', 'C[C>G]G': '7.75503247322284e-05', 'C[C>G]T': '9.56620779922714e-05', 'C[C>T]A': '0.0500324675691796', 'C[C>T]C': '0.0553359091315126', 'C[C>T]G': '0.0067143571477839', 'C[C>T]T': '0.0566367532883113', 'C[T>A]A': '0.000177114935194896', 'C[T>A]C': '0.000261169480711117', 'C[T>A]G': '0.000112072727354962', 'C[T>A]T': '0.000371240909363312', 'C[T>C]A': '0.000451292857474', 'C[T>C]C': '0.00025916818200835', 'C[T>C]G': '0.000104067532543894', 'C[T>C]T': '0.00182118181951814', 'C[T>G]A': '8.42546753864984e-05', 'C[T>G]C': '0.000172111688437978', 'C[T>G]G': '2.23144805358541e-16', 'C[T>G]T': '2.65172078116652e-05', 'G[C>A]A': '4.06263636661738e-05', 'G[C>A]C': '0.000138089610490936', 'G[C>A]G': '4.25275974338026e-05', 'G[C>A]T': '2.23144805358541e-16', 'G[C>G]A': '1.96127272871184e-05', 'G[C>G]C': '3.32215584659352e-05', 'G[C>G]G': '4.51292857474e-05', 'G[C>G]T': '4.20272727581108e-05', 'G[C>T]A': '8.49551299324669e-05', 'G[C>T]C': '0.011707597411188', 'G[C>T]G': '5.52358441963743e-06', 'G[C>T]T': '0.0161104545572758', 'G[T>A]A': '0.00040126038990482', 'G[T>A]C': '0.000224145454709925', 'G[T>A]G': '0.000290188311901242', 'G[T>A]

#     new_nt_seq = protein_nt_seq

#     for i in mutation_signatures.keys():
#         orig_tn = i[0] + i[2] + i[6]
#         new_tn = i[0] + i[4] + i[6]

#         prob = float(mutation_signatures[i])
        
#         if random.random() < (prob * 100):
#             new_nt_seq = new_nt_seq.replace(orig_tn, new_tn)


#     count = sum(1 for a, b in zip(protein_nt_seq, new_nt_seq) if a != b)

#     mutation_rate = count / len(protein_nt_seq)
#     return mutation_rate

# import os

# profile_directory = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/profiles'
# protein_directory = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/proteins'

# mut_rates = {}

# for root, dirs, files in os.walk(root_directory):
#     for profile in files:
#         print(profile)
#         # Construct the full file path
#         full_path = os.path.join(root, profile)
#         print(f"Processing file: {full_path}")
#         with open(full_path, 'r') as f:
#             content = f.read()
#         mut_rates[filename] = get_mutation_rate(protein, full_path)

# highest_key = max(mut_rates, key=mut_rates.get)

# # Get the highest value using the key
# highest_value = mut_rates[highest_key]

# print(highest_value, highest_key)

import os
import random

def get_mutation_rate(protein_nt_seq, file_path):
    mutation_signatures = {}
    
    with open(file_path, 'r') as f:
        # 1. Universally skip the first line (the header), no matter what it says
        next(f, None) 
        
        for line in f:
            parts = line.strip().split()
            # 2. Ensure the line isn't blank and has our data
            if len(parts) >= 2: 
                # Store the trinucleotide context (e.g., A[C>A]A) and its probability
                mutation_signatures[parts[0]] = float(parts[1])

    seq_list = list(protein_nt_seq.upper())
    mutations_count = 0

    # 3. Sliding window over the sequence to check every trinucleotide context
    for idx in range(len(seq_list) - 2):
        context = "".join(seq_list[idx:idx+3])

        # Check this context against ALL signature probabilities
        for mut_sig, prob in mutation_signatures.items():
            orig_tn = mut_sig[0] + mut_sig[2] + mut_sig[6] # e.g., A[C>A]A -> ACA
            
            if context == orig_tn:
                # Roll the dice using the extracted probability
                if random.random() < prob:
                    new_base = mut_sig[4] # Extract the mutated central base
                    seq_list[idx+1] = new_base # Apply mutation
                    mutations_count += 1
                    break # Once mutated, break out and move to the next position

    # Calculate final rate
    mutation_rate = mutations_count / len(protein_nt_seq) if len(protein_nt_seq) > 0 else 0
    return mutation_rate

# --- FILE STRUCTURE & ANALYSIS ---

profile_directory = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/profiles'
protein_directory = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/proteins'

# Your 4 target categories
categories = ['Tobacco', 'UV', 'Mismatch', 'Chemotherapy']
mut_rates = {}

print("Starting mutation simulation analysis...\n" + "-"*40)

for category in categories:
    cat_profile_dir = os.path.join(profile_directory, category)
    cat_protein_dir = os.path.join(protein_directory, category)

    if not os.path.exists(cat_profile_dir) or not os.path.exists(cat_protein_dir):
        print(f"Skipping '{category}': Subdirectories not found.")
        continue

    profiles = [f for f in os.listdir(cat_profile_dir) if not f.startswith('.')]
    proteins = [f for f in os.listdir(cat_protein_dir) if not f.startswith('.')]

    for profile in profiles:
        profile_path = os.path.join(cat_profile_dir, profile)

        for protein_file in proteins:
            protein_path = os.path.join(cat_protein_dir, protein_file)

            # Read protein sequence and remove FASTA headers
            with open(protein_path, 'r') as f:
                lines = f.readlines()
                protein_seq = "".join([line.strip() for line in lines if not line.startswith('>')])

            # Run the simulation for ALL signatures in the file
            rate = get_mutation_rate(protein_seq, profile_path)
            
            result_key = f"{category.upper()} | Profile: {profile} | Protein: {protein_file}"
            mut_rates[result_key] = rate
            print(f"Processed: {result_key} --> Rate: {rate:.6f}")

print("-" * 40)

if mut_rates:
    highest_key = max(mut_rates, key=mut_rates.get)
    highest_value = mut_rates[highest_key]
    print(f"\n🏆 HIGHEST MUTATION RATE FOUND:")
    print(f"{highest_key} \nRate: {highest_value:.6f} mutations/base")


