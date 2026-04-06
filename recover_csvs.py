"""
Recovery script — run this in the SAME Python session as your simulation
if the CSV saving failed at the end. It reconstructs all_results from the
printed output and saves the three CSV files without re-running anything.

Paste the printed leaderboard output into the PRINTED_OUTPUT string below,
or just call save_from_printed_output() with the path to a saved log file.

Alternatively, if your Python session is still alive and all_results is
still in memory, just call the three save functions directly:

    from protein_pair_analysis import (
        save_detailed_results_csv,
        save_cross_pair_summary_csv,
        save_wide_pivot_csv,
    )
    import os, numpy as np
    output_dir = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/results_recovered'
    os.makedirs(output_dir, exist_ok=True)
    save_detailed_results_csv(all_results, output_dir)
    save_cross_pair_summary_csv(cross_pair_avg, sig_type_map, null_cross_avg, output_dir)
    save_wide_pivot_csv(all_results, output_dir)
    print('Done!')
"""

import os
import csv
import re
import numpy as np

# ----------------------------------------------------------------
# OPTION A: If your Python session is still alive, use the comment
#           block above — it's the fastest path.
#
# OPTION B: If the session is closed, this script parses the printed
#           terminal output to reconstruct results and write CSVs.
#
# Point LOG_FILE at a saved copy of your terminal output,
# or paste the output directly into the string below.
# ----------------------------------------------------------------

LOG_FILE = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/prism_log.txt'

PRINTED_OUTPUT = """
PASTE YOUR TERMINAL OUTPUT HERE IF NOT USING LOG_FILE
"""

OUTPUT_DIR = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM/results_recovered'


def parse_log(text):
    """
    Parses printed PRISM output into all_results dict.
    Expects blocks like:
      PROTEIN PAIR: PAIR1-PDK
      ...
      --- Signature: SBS-MS (SBS8_PROFILE) [SBS] ---
      Final distance: 0.077510 | Movement toward P2: -0.050880
    """
    all_results = {}
    current_pair = None

    pair_re    = re.compile(r'PROTEIN PAIR:\s+(\S+)')
    sig_re     = re.compile(r'--- Signature:\s+(.+?)\s+\[(\w+)\]\s+---')
    move_re    = re.compile(r'Movement toward P2:\s+([-\d.]+)')
    null_re    = re.compile(r'Null model avg movement:\s+([-\d.]+)')

    lines = text.splitlines()
    current_sig = None
    current_type = None

    for line in lines:
        m = pair_re.search(line)
        if m:
            current_pair = m.group(1)
            all_results[current_pair] = {}
            continue

        m = sig_re.search(line)
        if m and current_pair:
            current_sig  = m.group(1).strip()
            current_type = m.group(2).strip()
            continue

        m = move_re.search(line)
        if m and current_pair and current_sig:
            net_change = float(m.group(1))
            all_results[current_pair][current_sig] = (net_change, current_type)
            current_sig = None
            continue

        m = null_re.search(line)
        if m and current_pair:
            null_val = float(m.group(1))
            all_results[current_pair]['~~NULL (uniform random)'] = (null_val, 'NULL')

    return all_results


def build_cross_pair(all_results):
    all_sig_names = set(
        sig for pr in all_results.values() for sig in pr
    )
    cross_pair_avg = {}
    sig_type_map   = {}
    for sig in all_sig_names:
        values = [all_results[p][sig][0] for p in all_results if sig in all_results[p]]
        for p in all_results:
            if sig in all_results[p]:
                sig_type_map[sig] = all_results[p][sig][1]
                break
        cross_pair_avg[sig] = (np.mean(values), np.std(values), len(values))
    null_cross_avg = cross_pair_avg.get('~~NULL (uniform random)', (0, 0, 0))[0]
    return cross_pair_avg, sig_type_map, null_cross_avg


def save_detailed(all_results, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, 'prism_detailed_results.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['protein_pair', 'signature', 'sig_type', 'net_movement_toward_p2'])
        for pair, pr in all_results.items():
            for sig, (val, stype) in pr.items():
                w.writerow([pair, sig, stype, f'{val:.6f}'])
    print(f'Saved: {path}')


def save_summary(cross_pair_avg, sig_type_map, null_cross_avg, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, 'prism_cross_pair_summary.csv')
    sorted_cross = sorted(cross_pair_avg.items(), key=lambda x: x[1][0], reverse=True)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'signature', 'sig_type', 'avg_movement_toward_p2',
                    'std', 'n_pairs', 'beats_null_model'])
        for rank, (sig, (avg, std, n)) in enumerate(sorted_cross, 1):
            beats = 'YES' if avg > null_cross_avg and sig != '~~NULL (uniform random)' else 'NO'
            w.writerow([rank, sig, sig_type_map.get(sig, '?'),
                        f'{avg:.6f}', f'{std:.6f}', n, beats])
    print(f'Saved: {path}')


def save_pivot(all_results, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, 'prism_pivot_table.csv')
    pairs = list(all_results.keys())
    all_sigs = sorted(set(sig for pr in all_results.values() for sig in pr))
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['signature'] + pairs + ['avg_across_pairs'])
        for sig in all_sigs:
            vals = [all_results[p][sig][0] if sig in all_results[p] else '' for p in pairs]
            nums = [v for v in vals if v != '']
            avg  = np.mean(nums) if nums else ''
            w.writerow([sig] + [f'{v:.6f}' if v != '' else '' for v in vals]
                       + [f'{avg:.6f}' if avg != '' else ''])
    print(f'Saved: {path}')


def main():
    text = open(LOG_FILE).read() if LOG_FILE else PRINTED_OUTPUT
    if not text.strip() or 'PASTE' in text:
        print("ERROR: No log text provided. Set LOG_FILE or paste output into PRINTED_OUTPUT.")
        return

    print("Parsing printed output...")
    all_results = parse_log(text)
    print(f"Found {len(all_results)} protein pairs:")
    for pair, pr in all_results.items():
        print(f"  {pair}: {len(pr)} signatures")

    cross_pair_avg, sig_type_map, null_cross_avg = build_cross_pair(all_results)

    save_detailed(all_results, OUTPUT_DIR)
    save_summary(cross_pair_avg, sig_type_map, null_cross_avg, OUTPUT_DIR)
    save_pivot(all_results, OUTPUT_DIR)
    print(f'\nAll CSVs saved to: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
