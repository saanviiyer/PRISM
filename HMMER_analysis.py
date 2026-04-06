"""
prism_hmmer.py  —  HMMER analysis integrated into the PRISM pipeline
======================================================================

DESIGN PHILOSOPHY
-----------------
Previous versions wrote empty CSVs because hmmscan / hmmsearch silently
skipped whenever Pfam-A.hmm or a sequence database was not present.
This rewrite fixes that by restructuring analyses into two tiers:

  TIER 1 — always runs, requires only HMMER binaries + your FASTA files
  -----------------------------------------------------------------------
  A. Pairwise profile scoring
       hmmbuild (P2 HMM) + hmmsearch (HMM vs P1) and vice versa
       hmmalign for structural coverage in both directions
       → bitscore, E-value, coverage — ALWAYS populated

  B. All-vs-all pairwise matrix
       One HMM per protein, every protein scored against every HMM
       → NxN bitscore matrix revealing sequence distance landscape

  TIER 2 — runs only when optional databases are present
  -----------------------------------------------------------------------
  C. Pfam domain annotation  (requires Pfam-A.hmm + hmmpress)
  D. Homolog search           (requires a local FASTA sequence database)

  INTEGRATION
  -----------
  Joins all HMMER results with PRISM cosine-distance CSVs when available.

PREREQUISITES
-------------
  brew install hmmer          # macOS
  apt install hmmer           # Linux
  conda install -c bioconda hmmer
  pip install biopython numpy

  # Optional Pfam:
  wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
  gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
  # Set PFAM_DB below.
"""

import os
import sys
import csv
import shutil
import subprocess
import textwrap
from datetime import datetime

import numpy as np


# =======================================================================
# CONFIGURATION
# =======================================================================

BASE_PATH = '/Users/saanviiyer/Downloads/CALTECH/RESEARCH/AITHYRA/PRISM'

# Tier 2 databases — set to '' to skip that analysis
PFAM_DB   = os.path.join(BASE_PATH, 'Pfam-A.hmm')
SEARCH_DB = ''
HMMBUILD  = 'hmmbuild'
HMMALIGN  = 'hmmalign'
HMMSCAN   = 'hmmscan'
HMMSEARCH = 'hmmsearch'

HMMER_CPUS     = 4
SCAN_EVALUE    = 1e-3
SEARCH_EVALUE  = 1e-5
TOP_N_HOMOLOGS = 20

PROTEIN_PAIRS = {
    'PAIR1-PDK': ('A0A010RJQ9.fa', 'A0A010S4M5.fa', 'PAIR1-PDK'),
    'PAIR2-Ras': ('A0A010R5C7.fa', 'A0A010SAT4.fa', 'PAIR2-Ras'),
    'PAIR3-p53': ('A0A087VQ67.fa', 'A0A087WT22.fa', 'PAIR3-p53'),
    'PAIR4-ABC': ('A0A010Z411.fa', 'A0A011NXC8.fa', 'PAIR4-ABC'),
    'PAIR5-HR':  ('A0A009T0V7.fa', 'A0A010QZB6.fa', 'PAIR5-HR'),
}

# Point at a PRISM detailed_results CSV to enable integration; else leave ''
PRISM_CSV = ''


# =======================================================================
# STARTUP CHECKS
# =======================================================================

def check_hmmer():
    missing = [b for b in [HMMBUILD, HMMALIGN, HMMSCAN, HMMSEARCH]
               if not shutil.which(b)]
    if missing:
        sys.exit(
            f"\nERROR: HMMER binaries not on PATH: {missing}\n"
            "Install: brew install hmmer | apt install hmmer | "
            "conda install -c bioconda hmmer\n"
        )
    print("HMMER binaries found. ✓")


def check_databases() -> dict:
    status = {
        'pfam':   os.path.exists(PFAM_DB),
        'search': os.path.exists(SEARCH_DB),
    }
    print(f"Pfam DB  : {'✓ ' + PFAM_DB   if status['pfam']   else '✗ NOT FOUND — domain annotation skipped'}")
    print(f"Search DB: {'✓ ' + SEARCH_DB  if status['search'] else '✗ NOT FOUND — homolog search skipped'}")

    if status['pfam']:
        missing_idx = [e for e in ['.h3i', '.h3m', '.h3f', '.h3p']
                       if not os.path.exists(PFAM_DB + e)]
        if missing_idx:
            print(f"  WARNING: Pfam index files missing: {missing_idx}")
            print(f"  Run:  hmmpress {PFAM_DB}")
            status['pfam'] = False
    return status


# =======================================================================
# FILE I/O
# =======================================================================

def read_fasta_seq(filepath: str) -> str:
    with open(filepath) as f:
        return ''.join(l.strip() for l in f if not l.startswith('>'))


def write_fasta(seq: str, seq_id: str, filepath: str):
    with open(filepath, 'w') as f:
        f.write(f'>{seq_id}\n')
        for chunk in textwrap.wrap(seq, 60):
            f.write(chunk + '\n')


def write_multi_fasta(entries: list, filepath: str):
    with open(filepath, 'w') as f:
        for seq_id, seq in entries:
            f.write(f'>{seq_id}\n')
            for chunk in textwrap.wrap(seq, 60):
                f.write(chunk + '\n')


# =======================================================================
# SUBPROCESS WRAPPER
# =======================================================================

def run_cmd(cmd: list, label: str):
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ✗ [{label}] exit {result.returncode}")
        for line in result.stderr.strip().splitlines():
            print(f"      {line}")
    return result.stdout, result.stderr, result.returncode


# =======================================================================
# PARSERS
# =======================================================================

def _sf(s):
    try:    return float(s)
    except: return None

def _si(s):
    try:    return int(s)
    except: return None


def parse_tblout(tbl_path: str, evalue_cutoff: float = 1.0,
                 top_n: int = 9999) -> list:
    """
    Parse hmmsearch --tblout.
    Columns (0-indexed): 0=target 1=t.acc 2=query 3=q.acc
      4=seq-Evalue 5=seq-score 6=seq-bias ... 18+=description
    """
    hits = []
    if not os.path.exists(tbl_path):
        print(f"      WARNING: tblout not found: {tbl_path}")
        return hits
    with open(tbl_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            try:
                ev = float(parts[4])
                sc = float(parts[5])
            except ValueError:
                continue
            if ev > evalue_cutoff:
                continue
            hits.append({
                'target':      parts[0],
                'query':       parts[2],
                'evalue':      ev,
                'bitscore':    sc,
                'bias':        _sf(parts[6]),
                'description': ' '.join(parts[18:]),
            })
            if len(hits) >= top_n:
                break
    return hits


def parse_domtblout(domtbl_path: str, evalue_cutoff: float = 1.0) -> list:
    """
    Parse hmmscan --domtblout.
    Columns: 0=domain 1=d.acc 2=tlen 3=query 4=q.acc 5=qlen
      6=seq-Eval 7=seq-score 8=seq-bias 9=dom# 10=#dom
      11=c-Eval 12=i-Eval(use this) 13=dom-score 14=dom-bias
      15=hmm-from 16=hmm-to 17=ali-from 18=ali-to
      19=env-from 20=env-to 21=meanPP 22+=description
    """
    hits = []
    if not os.path.exists(domtbl_path):
        print(f"      WARNING: domtblout not found: {domtbl_path}")
        return hits
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            try:
                i_ev  = float(parts[12])
                seq_ev = float(parts[6])
            except ValueError:
                continue
            if i_ev > evalue_cutoff:
                continue
            hits.append({
                'domain':      parts[0],
                'query':       parts[3],
                'seq_evalue':  seq_ev,
                'seq_score':   _sf(parts[7]),
                'dom_evalue':  i_ev,
                'dom_score':   _sf(parts[13]),
                'hmm_from':    _si(parts[15]),
                'hmm_to':      _si(parts[16]),
                'ali_from':    _si(parts[17]),
                'ali_to':      _si(parts[18]),
                'description': ' '.join(parts[22:]),
            })
    return hits


def parse_hmmbuild_stats(stdout: str) -> dict:
    """
    Parse hmmbuild stdout for model length and effective sequence count.
    hmmbuild prints a summary table like:
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # idx name                  nseq  alen  mlen eff_nseq re/pos
      # --- -------------------- ----- ----- ----- -------- ------
      #   1 PAIR1-PDK_P2             1   238   238     0.62  0.589
    The columns are fixed-width after the header. We parse the last
    data row (there is usually only one for single-sequence builds).
    Also handles older hmmbuild format that prints key: value lines.
    """
    stats = {}

    # Strategy 1: table format (hmmbuild 3.x default)
    # Header line contains 'idx', data lines start with '#   N'
    header_found = False
    for line in stdout.splitlines():
        s = line.strip()
        # Detect the column header line
        if s.startswith('# idx') or ('nseq' in s and 'alen' in s and 'mlen' in s):
            header_found = True
            continue
        # Data rows immediately follow the dashed separator after the header
        if header_found and s.startswith('#') and not s.startswith('# -'):
            # Strip leading '# ' and split
            parts = s.lstrip('#').split()
            # Expected: idx  name  nseq  alen  mlen  eff_nseq  re/pos
            if len(parts) >= 6:
                try:
                    stats['model_length'] = int(parts[4])    # mlen column
                    stats['neff']         = float(parts[5])  # eff_nseq column
                except (ValueError, IndexError):
                    pass

    # Strategy 2: key:value lines (some hmmbuild versions / verbose mode)
    if 'model_length' not in stats:
        for line in stdout.splitlines():
            s = line.strip()
            if 'Model length' in s or 'model length' in s:
                try:    stats['model_length'] = int(s.split()[-1])
                except: pass
            if 'Effective number' in s or 'eff_nseq' in s or 'Neff' in s:
                try:    stats['neff'] = float(s.split()[-1])
                except: pass

    return stats


def parse_hmmalign_coverage(sto_path: str) -> float:
    """Fraction of positions in the alignment that are residues (not gaps)."""
    aligned = gaps = 0
    if not os.path.exists(sto_path):
        return 0.0
    with open(sto_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('//') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            seq_str = parts[1]
            aligned += sum(1 for c in seq_str if c.isalpha())
            gaps    += seq_str.count('-')
    total = aligned + gaps
    return round(aligned / total, 4) if total > 0 else 0.0


# =======================================================================
# TIER 1A — PAIRWISE PROFILE SCORING (always runs)
# =======================================================================

def run_pairwise_scoring(pair_label: str, p1_path: str, p2_path: str,
                         work_dir: str) -> dict:
    """
    Build HMM from P2 → score P1 against it (and vice versa).
    hmmalign gives structural coverage.
    All fields are always written; None means the run failed.
    """
    print(f"\n  [A] Pairwise scoring: {pair_label}")
    result = {
        'pair':              pair_label,
        'model_length_p2':   None, 'model_length_p1': None,
        'neff_p2':           None, 'neff_p1':          None,
        'bitscore_p1_vs_p2': None, 'evalue_p1_vs_p2':  None,
        'bitscore_p2_vs_p1': None, 'evalue_p2_vs_p1':  None,
        'p1_coverage_on_p2': None, 'p2_coverage_on_p1':None,
        'error':             None,
    }

    p1_seq = read_fasta_seq(p1_path)
    p2_seq = read_fasta_seq(p2_path)
    p1_fa  = os.path.join(work_dir, f'{pair_label}_p1.fa')
    p2_fa  = os.path.join(work_dir, f'{pair_label}_p2.fa')
    hmm_p1 = os.path.join(work_dir, f'{pair_label}_p1.hmm')
    hmm_p2 = os.path.join(work_dir, f'{pair_label}_p2.hmm')
    write_fasta(p1_seq, f'{pair_label}_P1', p1_fa)
    write_fasta(p2_seq, f'{pair_label}_P2', p2_fa)

    errors = []

    # Build P2 HMM
    stdout, _, rc = run_cmd(
        [HMMBUILD, '--cpu', str(HMMER_CPUS), hmm_p2, p2_fa], 'hmmbuild P2')
    if rc != 0:
        errors.append('hmmbuild-P2')
    else:
        st = parse_hmmbuild_stats(stdout)
        result['model_length_p2'] = st.get('model_length')
        result['neff_p2']         = st.get('neff')

    # Build P1 HMM
    stdout, _, rc = run_cmd(
        [HMMBUILD, '--cpu', str(HMMER_CPUS), hmm_p1, p1_fa], 'hmmbuild P1')
    if rc != 0:
        errors.append('hmmbuild-P1')
    else:
        st = parse_hmmbuild_stats(stdout)
        result['model_length_p1'] = st.get('model_length')
        result['neff_p1']         = st.get('neff')

    # Score P1 against P2-profile
    if 'hmmbuild-P2' not in errors:
        tbl = os.path.join(work_dir, f'{pair_label}_p1_vs_p2_tbl.txt')
        out = os.path.join(work_dir, f'{pair_label}_p1_vs_p2.txt')
        run_cmd([HMMSEARCH, '--cpu', str(HMMER_CPUS),
                 '--tblout', tbl, '-E', '1e6', '--max', '--noali',
                 '-o', out, hmm_p2, p1_fa],
                'hmmsearch P1 vs P2-profile')
        hits = parse_tblout(tbl, evalue_cutoff=1e6)
        if hits:
            result['bitscore_p1_vs_p2'] = hits[0]['bitscore']
            result['evalue_p1_vs_p2']   = hits[0]['evalue']
        else:
            print("      No hit returned — sequences may be highly divergent")

        # hmmalign coverage P1 on P2-profile
        sto = os.path.join(work_dir, f'{pair_label}_p1_on_p2.sto')
        run_cmd([HMMALIGN, '--trim', '--outformat', 'stockholm',
                 '-o', sto, hmm_p2, p1_fa],
                'hmmalign P1 on P2-profile')
        result['p1_coverage_on_p2'] = parse_hmmalign_coverage(sto)

    # Score P2 against P1-profile (reciprocal)
    if 'hmmbuild-P1' not in errors:
        tbl = os.path.join(work_dir, f'{pair_label}_p2_vs_p1_tbl.txt')
        out = os.path.join(work_dir, f'{pair_label}_p2_vs_p1.txt')
        run_cmd([HMMSEARCH, '--cpu', str(HMMER_CPUS),
                 '--tblout', tbl, '-E', '1e6', '--max', '--noali',
                 '-o', out, hmm_p1, p2_fa],
                'hmmsearch P2 vs P1-profile')
        hits = parse_tblout(tbl, evalue_cutoff=1e6)
        if hits:
            result['bitscore_p2_vs_p1'] = hits[0]['bitscore']
            result['evalue_p2_vs_p1']   = hits[0]['evalue']

        sto = os.path.join(work_dir, f'{pair_label}_p2_on_p1.sto')
        run_cmd([HMMALIGN, '--trim', '--outformat', 'stockholm',
                 '-o', sto, hmm_p1, p2_fa],
                'hmmalign P2 on P1-profile')
        result['p2_coverage_on_p1'] = parse_hmmalign_coverage(sto)

    if errors:
        result['error'] = '; '.join(errors)

    print(f"    P1→P2  bitscore={result['bitscore_p1_vs_p2']}  "
          f"E={result['evalue_p1_vs_p2']}  cov={result['p1_coverage_on_p2']}")
    print(f"    P2→P1  bitscore={result['bitscore_p2_vs_p1']}  "
          f"E={result['evalue_p2_vs_p1']}  cov={result['p2_coverage_on_p1']}")
    return result


# =======================================================================
# TIER 1B — ALL-VS-ALL MATRIX (always runs)
# =======================================================================

def run_all_vs_all(all_seqs: dict, work_dir: str):
    """
    Build one HMM per sequence, score every sequence against every HMM.
    Returns (matrix_dict, sorted_labels).
    matrix_dict: { (query_label, target_hmm_label): {'bitscore', 'evalue'} }
    """
    print(f"\n  [B] All-vs-all matrix ({len(all_seqs)} sequences)")
    labels  = sorted(all_seqs.keys())
    hmm_dir = os.path.join(work_dir, 'hmms')
    fa_dir  = os.path.join(work_dir, 'fas')
    os.makedirs(hmm_dir, exist_ok=True)
    os.makedirs(fa_dir,  exist_ok=True)

    # One FASTA + HMM per sequence
    hmms = {}
    for label, seq in all_seqs.items():
        fa  = os.path.join(fa_dir,  f'{label}.fa')
        hmm = os.path.join(hmm_dir, f'{label}.hmm')
        write_fasta(seq, label, fa)
        _, _, rc = run_cmd(
            [HMMBUILD, '--cpu', str(HMMER_CPUS), hmm, fa], f'hmmbuild {label}')
        if rc == 0:
            hmms[label] = hmm

    # One multi-FASTA of all queries (reused for every hmmsearch call)
    all_fa = os.path.join(fa_dir, '_all.fa')
    write_multi_fasta([(lbl, all_seqs[lbl]) for lbl in labels], all_fa)

    matrix = {}
    for target_label, hmm_path in hmms.items():
        tbl = os.path.join(work_dir, f'ava_{target_label}_tbl.txt')
        out = os.path.join(work_dir, f'ava_{target_label}.txt')
        run_cmd([HMMSEARCH, '--cpu', str(HMMER_CPUS),
                 '--tblout', tbl,
                 '-E', '1e6',    # very permissive — catch even remote homologs
                 '--max',        # disable all heuristic filters (no bias filter,
                                 # no composition correction) so divergent sequences
                                 # are not silently dropped
                 '--noali',
                 '-o', out, hmm_path, all_fa],
                f'ava vs {target_label}')
        for hit in parse_tblout(tbl, evalue_cutoff=1e6):
            matrix[(hit['target'], target_label)] = {
                'bitscore': hit['bitscore'],
                'evalue':   hit['evalue'],
            }

    filled = sum(1 for v in matrix.values() if v)
    print(f"    {filled} / {len(labels)**2} cells populated")
    return matrix, labels


# =======================================================================
# TIER 2A — DOMAIN ANNOTATION (optional, requires Pfam)
# =======================================================================

def run_domain_annotation(pair_label: str, p1_path: str, p2_path: str,
                           work_dir: str) -> dict:
    """
    Pfam domain annotation — two modes:
      1. Local hmmscan  (fast, requires Pfam-A.hmm + hmmpress)
      2. EBI HMMER API  (no local DB needed, ~5-10s per sequence via HTTPS)
    Mode 2 is used automatically when PFAM_DB is not present.
    """
    result = {
        'pair': pair_label,
        'p1_domains': [], 'p2_domains': [],
        'p1_domain_hits': [], 'p2_domain_hits': [],
        'shared_domains': [], 'domain_overlap_fraction': 0.0,
        'skipped': False,
        'source': '',
    }

    if os.path.exists(PFAM_DB):
        result['source'] = 'local_hmmscan'
        print(f"\n  [C] Domain annotation (local Pfam): {pair_label}")
        for label, path, key in [('P1', p1_path, 'p1'), ('P2', p2_path, 'p2')]:
            seq    = read_fasta_seq(path)
            fa     = os.path.join(work_dir, f'{pair_label}_{key}_scan.fa')
            domtbl = os.path.join(work_dir, f'{pair_label}_{key}_domtbl.txt')
            out    = os.path.join(work_dir, f'{pair_label}_{key}_scan.txt')
            write_fasta(seq, f'{pair_label}_{label}', fa)
            run_cmd([HMMSCAN, '--cpu', str(HMMER_CPUS),
                     '--domtblout', domtbl,
                     '-E', str(SCAN_EVALUE), '--domE', str(SCAN_EVALUE),
                     '--noali', '-o', out, PFAM_DB, fa],
                    f'hmmscan {label}')
            hits    = parse_domtblout(domtbl, SCAN_EVALUE)
            domains = sorted({h['domain'] for h in hits})
            result[f'{key}_domains']     = domains
            result[f'{key}_domain_hits'] = hits
            print(f"    {label}: {len(domains)} domains  "
                  + str(domains[:4]) + (' ...' if len(domains) > 4 else ''))
    else:
        # ---- EBI HMMER API fallback ----
        try:
            import urllib.request, urllib.parse, time, json as _json
        except ImportError:
            print(f"\n  [C] Domain annotation SKIPPED — no Pfam DB and urllib unavailable")
            result['skipped'] = True
            return result

        result['source'] = 'ebi_api'
        print(f"\n  [C] Domain annotation via EBI HMMER API: {pair_label}")
        EBI_URL = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'

        def _ebi_pfam_scan(seq: str, label: str) -> list:
            """POST sequence to EBI hmmscan, poll until done, return hit list."""
            data = urllib.parse.urlencode({
                'seqdb': 'pfam',
                'seq':   seq,
            }).encode()
            req = urllib.request.Request(
                EBI_URL, data=data,
                headers={'Accept': 'application/json',
                         'Content-Type': 'application/x-www-form-urlencoded'})
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    result_url = resp.geturl()   # EBI redirects to result page
                    raw = resp.read().decode()
            except Exception as e:
                print(f"      EBI API error for {label}: {e}")
                return []

            try:
                payload = _json.loads(raw)
            except Exception:
                print(f"      EBI API: could not parse JSON response for {label}")
                return []

            hits = []
            for hit in payload.get('results', {}).get('hits', []):
                for dom in hit.get('domains', []):
                    i_ev = dom.get('ievalue', 1.0)
                    if float(i_ev) > SCAN_EVALUE:
                        continue
                    hits.append({
                        'domain':      hit.get('name', ''),
                        'query':       label,
                        'seq_evalue':  hit.get('evalue', None),
                        'seq_score':   hit.get('score', None),
                        'dom_evalue':  i_ev,
                        'dom_score':   dom.get('bitscore', None),
                        'hmm_from':    dom.get('hmm_from', None),
                        'hmm_to':      dom.get('hmm_to', None),
                        'ali_from':    dom.get('ali_from', None),
                        'ali_to':      dom.get('ali_to', None),
                        'description': hit.get('desc', ''),
                    })
            return hits

        for label, path, key in [('P1', p1_path, 'p1'), ('P2', p2_path, 'p2')]:
            seq  = read_fasta_seq(path)
            print(f"    Querying EBI for {label} ({len(seq)} aa) ...")
            hits = _ebi_pfam_scan(seq, f'{pair_label}_{label}')
            domains = sorted({h['domain'] for h in hits})
            result[f'{key}_domains']     = domains
            result[f'{key}_domain_hits'] = hits
            print(f"    {label}: {len(domains)} domains  "
                  + str(domains[:4]) + (' ...' if len(domains) > 4 else ''))
            import time; time.sleep(1)   # be polite to EBI

    # Shared domain Jaccard
    p1s = set(result['p1_domains'])
    p2s = set(result['p2_domains'])
    shared = p1s & p2s
    union  = p1s | p2s
    result['shared_domains']          = sorted(shared)
    result['domain_overlap_fraction'] = round(len(shared) / len(union), 4) if union else 0.0
    print(f"    Shared: {sorted(shared)}  Jaccard={result['domain_overlap_fraction']:.3f}")
    return result


# =======================================================================
# TIER 2B — HOMOLOG SEARCH (optional, requires SEARCH_DB)
# =======================================================================

def run_homolog_search(pair_label: str, p1_path: str, work_dir: str) -> dict:
    result = {'pair': pair_label, 'hits': [], 'n_hits': 0,
              'skipped': not os.path.exists(SEARCH_DB)}
    if result['skipped']:
        return result

    print(f"\n  [D] Homolog search: {pair_label}")
    seq = read_fasta_seq(p1_path)
    fa  = os.path.join(work_dir, f'{pair_label}_p1_search.fa')
    hmm = os.path.join(work_dir, f'{pair_label}_p1_search.hmm')
    tbl = os.path.join(work_dir, f'{pair_label}_search_tbl.txt')
    out = os.path.join(work_dir, f'{pair_label}_search.txt')
    write_fasta(seq, f'{pair_label}_P1', fa)
    run_cmd([HMMBUILD, '--cpu', str(HMMER_CPUS), hmm, fa], 'hmmbuild P1')
    run_cmd([HMMSEARCH, '--cpu', str(HMMER_CPUS),
             '--tblout', tbl, '-E', str(SEARCH_EVALUE),
             '--noali', '-o', out, hmm, SEARCH_DB],
            'hmmsearch vs DB')
    hits = parse_tblout(tbl, SEARCH_EVALUE, TOP_N_HOMOLOGS)
    result['hits']   = hits
    result['n_hits'] = len(hits)
    if hits:
        b = hits[0]
        print(f"    Top: {b['target']}  E={b['evalue']:.2e}  score={b['bitscore']:.1f}")
    else:
        print("    No significant homologs found.")
    return result


# =======================================================================
# INTEGRATION
# =======================================================================

def integrate_with_prism(scoring: list, domains: list, prism_csv: str) -> list:
    if not prism_csv or not os.path.exists(prism_csv):
        return []
    s_map = {r['pair']: r for r in scoring}
    d_map = {r['pair']: r for r in domains if not r.get('skipped')}
    rows = []
    with open(prism_csv, newline='') as f:
        for row in csv.DictReader(f):
            pair = row.get('protein_pair', '')
            s    = s_map.get(pair, {})
            d    = d_map.get(pair, {})
            rows.append({
                'protein_pair':            pair,
                'signature':               row.get('signature', ''),
                'sig_type':                row.get('sig_type', ''),
                'net_movement_toward_p2':  row.get('net_movement_toward_p2', ''),
                'hmmer_bitscore_p1_vs_p2': s.get('bitscore_p1_vs_p2', ''),
                'hmmer_evalue_p1_vs_p2':   s.get('evalue_p1_vs_p2', ''),
                'hmmer_bitscore_p2_vs_p1': s.get('bitscore_p2_vs_p1', ''),
                'hmmer_p1_coverage':       s.get('p1_coverage_on_p2', ''),
                'hmmer_p2_coverage':       s.get('p2_coverage_on_p1', ''),
                'pfam_jaccard':            d.get('domain_overlap_fraction', ''),
                'pfam_shared_domains':     '|'.join(d.get('shared_domains', [])),
            })
    print(f"    Integrated {len(rows)} rows")
    return rows


# =======================================================================
# CSV SAVERS
# =======================================================================

def save_scoring_csv(results: list, output_dir: str):
    path = os.path.join(output_dir, 'hmmer_pairwise_scoring.csv')
    fields = ['pair',
              'model_length_p2', 'model_length_p1', 'neff_p2', 'neff_p1',
              'bitscore_p1_vs_p2', 'evalue_p1_vs_p2',
              'bitscore_p2_vs_p1', 'evalue_p2_vs_p1',
              'p1_coverage_on_p2', 'p2_coverage_on_p1', 'error']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(results)
    print(f"  → {path}")


def save_ava_csv(matrix: dict, labels: list, output_dir: str):
    path = os.path.join(output_dir, 'hmmer_allvsall_bitscore.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['query \\ target'] + labels)
        for q in labels:
            row = [q]
            for t in labels:
                cell = matrix.get((q, t))
                row.append(f"{cell['bitscore']:.1f}" if cell else '')
            w.writerow(row)
    print(f"  → {path}")

    # Also save E-value matrix
    path_ev = os.path.join(output_dir, 'hmmer_allvsall_evalue.csv')
    with open(path_ev, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['query \\ target'] + labels)
        for q in labels:
            row = [q]
            for t in labels:
                cell = matrix.get((q, t))
                row.append(f"{cell['evalue']:.2e}" if cell else '')
            w.writerow(row)
    print(f"  → {path_ev}")


def save_domain_csvs(results: list, output_dir: str):
    # Per-domain hits
    path = os.path.join(output_dir, 'hmmer_domain_hits.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['pair', 'sequence', 'pfam_domain', 'dom_evalue',
                    'dom_score', 'ali_from', 'ali_to', 'description'])
        for r in results:
            if r.get('skipped'):
                continue
            for key, label in [('p1_domain_hits', 'P1'), ('p2_domain_hits', 'P2')]:
                for h in r.get(key, []):
                    w.writerow([r['pair'], label, h['domain'], h['dom_evalue'],
                                h['dom_score'], h['ali_from'], h['ali_to'],
                                h['description']])
    print(f"  → {path}")

    # Pair summary
    path2 = os.path.join(output_dir, 'hmmer_domain_summary.csv')
    with open(path2, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['pair', 'source', 'n_p1_domains', 'n_p2_domains', 'n_shared',
                    'domain_overlap_jaccard', 'p1_domains', 'p2_domains',
                    'shared_domains'])
        for r in results:
            if r.get('skipped'):
                w.writerow([r['pair'], 'SKIPPED'] + [''] * 7)
            else:
                w.writerow([r['pair'],
                            r.get('source', ''),
                            len(r['p1_domains']), len(r['p2_domains']),
                            len(r['shared_domains']),
                            r['domain_overlap_fraction'],
                            '|'.join(r['p1_domains']),
                            '|'.join(r['p2_domains']),
                            '|'.join(r['shared_domains'])])
    print(f"  → {path2}")


def save_homolog_csv(results: list, output_dir: str):
    path = os.path.join(output_dir, 'hmmer_homologs.csv')
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['pair', 'rank', 'target', 'evalue', 'bitscore', 'description'])
        for r in results:
            if r.get('skipped'):
                w.writerow([r['pair'], 'SKIPPED', '', '', '', ''])
            elif not r['hits']:
                w.writerow([r['pair'], 'no_hits', '', '', '', ''])
            else:
                for rank, h in enumerate(r['hits'], 1):
                    w.writerow([r['pair'], rank, h['target'],
                                h['evalue'], h['bitscore'], h['description']])
    print(f"  → {path}")


def save_integrated_csv(rows: list, output_dir: str):
    if not rows:
        return
    path = os.path.join(output_dir, 'hmmer_prism_integrated.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  → {path}")


# =======================================================================
# MAIN
# =======================================================================

def main():
    check_hmmer()
    check_databases()

    timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(BASE_PATH, f'hmmer_results_{timestamp}')
    work_dir   = os.path.join(output_dir, '_tmp')
    os.makedirs(work_dir, exist_ok=True)
    print(f"\nOutput : {output_dir}\n")

    all_seqs   = {}
    pair_paths = {}

    for pair_label, (p1_file, p2_file, subfolder) in PROTEIN_PAIRS.items():
        p1 = os.path.join(BASE_PATH, 'protein_pairs', subfolder, p1_file)
        p2 = os.path.join(BASE_PATH, 'protein_pairs', subfolder, p2_file)
        if not os.path.exists(p1) or not os.path.exists(p2):
            print(f"  WARNING: missing FASTA for {pair_label} — skipping")
            continue
        pair_paths[pair_label] = (p1, p2)
        all_seqs[f'{pair_label}_P1'] = read_fasta_seq(p1)
        all_seqs[f'{pair_label}_P2'] = read_fasta_seq(p2)

    scoring_results = []
    domain_results  = []
    homolog_results = []

    for pair_label, (p1, p2) in pair_paths.items():
        print(f"\n{'='*70}\nPAIR: {pair_label}\n{'='*70}")
        pw = os.path.join(work_dir, pair_label)
        os.makedirs(pw, exist_ok=True)

        scoring_results.append(run_pairwise_scoring(pair_label, p1, p2, pw))
        domain_results.append(run_domain_annotation(pair_label, p1, p2, pw))
        homolog_results.append(run_homolog_search(pair_label, p1, pw))

    ava_work = os.path.join(work_dir, '_ava')
    os.makedirs(ava_work, exist_ok=True)
    matrix, labels = run_all_vs_all(all_seqs, ava_work)

    print(f"\n{'='*70}\nSAVING\n{'='*70}")
    save_scoring_csv(scoring_results, output_dir)
    save_ava_csv(matrix, labels, output_dir)
    save_domain_csvs(domain_results, output_dir)
    save_homolog_csv(homolog_results, output_dir)
    integrated = integrate_with_prism(scoring_results, domain_results, PRISM_CSV)
    if integrated:
        save_integrated_csv(integrated, output_dir)
    else:
        print("  Integration skipped — set PRISM_CSV at top of file to enable")

    # Console summary
    print(f"\n{'='*70}\nSUMMARY\n{'='*70}")
    print(f"{'Pair':<15} {'BS P1→P2':>10} {'Cov P1':>8} "
          f"{'BS P2→P1':>10} {'Cov P2':>8} {'Domains':>8} {'Homologs':>9}")
    print('-' * 75)
    s_map = {r['pair']: r for r in scoring_results}
    d_map = {r['pair']: r for r in domain_results}
    h_map = {r['pair']: r for r in homolog_results}
    for pair in pair_paths:
        s = s_map.get(pair, {})
        d = d_map.get(pair, {})
        h = h_map.get(pair, {})
        print(
            f"{pair:<15}"
            f" {str(s.get('bitscore_p1_vs_p2','N/A')):>10}"
            f" {str(s.get('p1_coverage_on_p2','N/A')):>8}"
            f" {str(s.get('bitscore_p2_vs_p1','N/A')):>10}"
            f" {str(s.get('p2_coverage_on_p1','N/A')):>8}"
            f" {'SKIP' if d.get('skipped') else str(len(d.get('shared_domains',[]))):>8}"
            f" {'SKIP' if h.get('skipped') else str(h.get('n_hits',0)):>9}"
        )

    print(f"\nAll files written to: {output_dir}")
    print("  hmmer_pairwise_scoring.csv   — always populated")
    print("  hmmer_allvsall_bitscore.csv  — always populated")
    print("  hmmer_allvsall_evalue.csv    — always populated")
    print("  hmmer_domain_hits.csv        — requires Pfam-A.hmm")
    print("  hmmer_domain_summary.csv     — requires Pfam-A.hmm")
    print("  hmmer_homologs.csv           — requires SEARCH_DB")


if __name__ == '__main__':
    main()