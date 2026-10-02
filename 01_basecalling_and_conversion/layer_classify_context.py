#!/usr/bin/env python3
"""Keep reference cytosines and classify each as CpG or CpH, strand-aware, from the reference.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, gzip, subprocess, pysam
from collections import Counter

pileup, regions, ref, prefix = sys.argv[1:5]
fa = pysam.FastaFile(ref)

proc = subprocess.Popen(['tabix', '-R', regions, pileup],
                        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

out_g = gzip.open(prefix + '.cpg.tsv.gz', 'wt', compresslevel=4)
out_h = gzip.open(prefix + '.cph.tsv.gz', 'wt', compresslevel=4)
for o in (out_g, out_h):
    o.write("chrom\tstart\tend\tstrand\tNvalid\tNmod\tcontext\n")

cur, seq, L = None, None, 0
stats = Counter()
COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

for line in proc.stdout:
    f = line.rstrip('\n').split('\t')
    if f[3] != 'm':
        stats['not_m'] += 1
        continue
    c = f[0]
    if c != cur:
        cur = c
        seq = fa.fetch(c).upper()
        L = len(seq)
    p = int(f[1]); st = f[5]
    if p <= 0 or p >= L - 1:
        stats['edge'] += 1
        continue
    b = seq[p]
    if st == '+':
        if b != 'C':
            stats['not_ref_C'] += 1
            continue
        nxt = seq[p + 1]
        ctx = 'CpG' if nxt == 'G' else 'Cp' + nxt
    else:
        if b != 'G':
            stats['not_ref_C'] += 1
            continue
        prv = seq[p - 1]

        ctx = 'CpG' if prv == 'C' else 'Cp' + COMP.get(prv, 'N')
    if ctx.endswith('N'):
        stats['N_context'] += 1
        continue
    row = f"{c}\t{f[1]}\t{f[2]}\t{st}\t{f[9]}\t{f[11]}\t{ctx}\n"
    if ctx == 'CpG':
        out_g.write(row); stats['CpG'] += 1
    else:
        out_h.write(row); stats['CpH'] += 1
        stats[ctx] += 1

out_g.close(); out_h.close(); proc.stdout.close(); proc.wait()

tot = stats['CpG'] + stats['CpH']
print(f"kept  : {tot:,}   CpG {stats['CpG']:,} ({100*stats['CpG']/tot:.2f}%)  "
      f"CpH {stats['CpH']:,} ({100*stats['CpH']/tot:.2f}%)")
print(f"        CpA {stats['CpA']:,}  CpT {stats['CpT']:,}  CpC {stats['CpC']:,}")
print(f"dropped: not-ref-C {stats['not_ref_C']:,}  "
      f"({100*stats['not_ref_C']/(tot+stats['not_ref_C']):.2f}% of C rows)  "
      f"edge {stats['edge']:,}  N-context {stats['N_context']:,}  non-m {stats['not_m']:,}")
