#!/usr/bin/env python3
"""Per-read methylation (beta) from read-level CpG calls.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, os, gzip
from collections import defaultdict

P = '/data/8oxo_project/AdMSC_DMR/readlevel_calls'
RES = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-11/results'
os.makedirs(RES, exist_ok=True)

NBINS = 200
THRESHOLDS = [0.70, 0.80, 0.90, 0.98]
DEEP = 50

ARMS = ['NT', 'PA']

def is_cpg(kmer, ref_strand):
    """CpG dyad, on either strand, from the reference-forward 5-mer centred on the call."""
    k = kmer.lower()
    if len(k) < 5:
        return False
    if ref_strand == '-':
        return k[2] == 'g' and k[1] == 'c'
    return k[2] == 'c' and k[3] == 'g'

def pass1_depth(arm):
    """count reads per (chrom, pos, strand) so the deep stratum can be defined"""
    depth = defaultdict(int)
    n = 0
    with open(f'{P}/{arm}.calls.tsv', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 16:
                continue
            if f[17] != 'C':
                continue
            if not is_cpg(f[15], f[5]):
                continue
            depth[(f[3], f[2], f[5])] += 1
            n += 1
    return depth, n

def pass2(arm, deep_sites):
    hist = {'CpG': [0] * (NBINS + 1), 'CpH': [0] * (NBINS + 1),
            'CpG_deep': [0] * (NBINS + 1)}

    reads_seen = set()
    reopened = 0
    cur = None
    percall = {t: [0, 0] for t in THRESHOLDS}
    beta_hist = {t: [0] * 21 for t in THRESHOLDS}
    calls_per_read = {t: defaultdict(int) for t in THRESHOLDS}
    nreads = {t: 0 for t in THRESHOLDS}
    pooled = {t: [0, 0] for t in THRESHOLDS}
    total = 0

    def flush():
        for t in THRESHOLDS:
            v, m = percall[t]
            if v >= 1:
                nreads[t] += 1
                calls_per_read[t][min(v, 60)] += 1
                beta_hist[t][min(20, int(round(20 * m / v)))] += 1
            percall[t][0] = percall[t][1] = 0

    with open(f'{P}/{arm}.calls.tsv', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 16:
                continue
            rid = f[0]
            if rid != cur:
                if cur is not None:
                    flush()
                    if rid in reads_seen:
                        reopened += 1
                    reads_seen.add(rid)
                cur = rid
            if f[17] != 'C':
                continue
            try:
                q = float(f[12])
            except ValueError:
                continue
            cpg = is_cpg(f[15], f[5])

            b = min(NBINS, int(q * NBINS))
            hist['CpG' if cpg else 'CpH'][b] += 1
            if cpg and (f[3], f[2], f[5]) in deep_sites:
                hist['CpG_deep'][b] += 1
            if not cpg:
                continue
            total += 1
            mod = (f[13] == 'm')
            win = q
            for t in THRESHOLDS:
                if win >= t:
                    percall[t][0] += 1
                    pooled[t][0] += 1
                    if mod:
                        percall[t][1] += 1
                        pooled[t][1] += 1
    if cur is not None:
        flush()
    return dict(hist=hist, beta=beta_hist, cpr=calls_per_read, nreads=nreads,
                pooled=pooled, total=total, reopened=reopened, nread_ids=len(reads_seen))

out = {}
for arm in ARMS:
    print(f'[{arm}] pass 1 — per-position depth', flush=True)
    depth, ncpg = pass1_depth(arm)
    deep = {k for k, v in depth.items() if v >= DEEP}
    print(f'[{arm}]   {ncpg:,} CpG calls at {len(depth):,} positions; '
          f'{len(deep):,} positions at depth >= {DEEP}', flush=True)
    print(f'[{arm}] pass 2 — distributions', flush=True)
    out[arm] = pass2(arm, deep)
    r = out[arm]
    print(f'[{arm}]   {r["nread_ids"]:,} reads, {r["reopened"]:,} read ids reappeared '
          f'(0 means the file is grouped by read)', flush=True)
    for t in THRESHOLDS:
        v, m = r['pooled'][t]
        print(f'[{arm}]   thr {t:.2f}: {v:,} valid CpG calls, '
              f'{100*m/v if v else 0:.4f} % modified, {r["nreads"][t]:,} reads', flush=True)

with open(f'{RES}/level0a_prob_histogram.tsv', 'w') as o:
    o.write("arm\tcontext\tbin_lo\tbin_hi\tn_calls\n")
    for arm in ARMS:
        for ctx, h in out[arm]['hist'].items():
            for i, c in enumerate(h):
                if c:
                    o.write(f"{arm}\t{ctx}\t{i/NBINS:.4f}\t{(i+1)/NBINS:.4f}\t{c}\n")

with open(f'{RES}/level0a_read_beta.tsv', 'w') as o:
    o.write("arm\tthreshold\tbeta_bin_lo\tbeta_bin_hi\tn_reads\n")
    for arm in ARMS:
        for t in THRESHOLDS:
            for i, c in enumerate(out[arm]['beta'][t]):
                o.write(f"{arm}\t{t}\t{i/20:.2f}\t{(i+1)/20:.2f}\t{c}\n")

with open(f'{RES}/level0a_calls_per_read.tsv', 'w') as o:
    o.write("arm\tthreshold\tcalls_in_read\tn_reads\n")
    for arm in ARMS:
        for t in THRESHOLDS:
            for k in sorted(out[arm]['cpr'][t]):
                o.write(f"{arm}\t{t}\t{k}\t{out[arm]['cpr'][t][k]}\n")

with open(f'{RES}/level0a_pooled_by_threshold.tsv', 'w') as o:
    o.write("threshold\tNT_valid\tNT_mod\tNT_pct\tPA_valid\tPA_mod\tPA_pct\tratio\n")
    for t in THRESHOLDS:
        nv, nm = out['NT']['pooled'][t]
        pv, pm = out['PA']['pooled'][t]
        r = (pm / pv) / (nm / nv) if nv and pv and nm else float('nan')
        o.write(f"{t}\t{nv}\t{nm}\t{100*nm/nv:.5f}\t{pv}\t{pm}\t{100*pm/pv:.5f}\t{r:.5f}\n")
        print(f"  READ-SPACE thr {t:.2f}:  NT {100*nm/nv:.4f} %   PA {100*pm/pv:.4f} %   "
              f"ratio {r:.5f}", flush=True)

print(f"\nwrote {RES}/level0a_*.tsv")
