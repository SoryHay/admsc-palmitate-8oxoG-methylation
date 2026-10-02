#!/usr/bin/env python3
"""Quality-decile split and aggregation (Table S4).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, csv, collections
D, out = sys.argv[1], sys.argv[2]
L = ['# CHECK — 5mC and 5hmC call fractions by decile of read mean quality (qs), within each arm\n',
     'Command: `p1_check_qdecile.sh` (modkit extract calls --cpg --filter-threshold C:0.7 --pass-only --mapped-only) → `p1_check_qdecile.py`. '
     'Reads only whether the m/h assignment tracks read quality inside one library; says nothing about whether the arm-to-arm difference is biological.\n']
for arm in ('NT', 'PA'):
    qs = {}
    for line in open(f'{D}/{arm}_qs.tsv'):
        r, q = line.rstrip('\n').split('\t'); qs[r] = float(q)
    per_read = collections.defaultdict(lambda: [0, 0, 0])
    with open(f'{D}/{arm}_calls_cpg.tsv') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for row in rd:
            rec = per_read[row['read_id']]; rec[0] += 1
            c = row['call_code']
            if c == 'm': rec[1] += 1
            elif c == 'h': rec[2] += 1
    reads = [(qs[r], v) for r, v in per_read.items() if r in qs]
    reads.sort(key=lambda t: t[0]); n = len(reads)
    L.append(f'\n## {arm} — {n:,} reads with CpG calls\n')
    L.append('| decile | Q range | reads | valid calls | 5mC % | 5hmC % | 5mC+5hmC % |\n|---|---|---|---|---|---|---|')
    tot = [0, 0, 0]
    for d in range(10):
        chunk = reads[d * n // 10:(d + 1) * n // 10]
        v = sum(c[1][0] for c in chunk); m = sum(c[1][1] for c in chunk); h = sum(c[1][2] for c in chunk)
        tot[0] += v; tot[1] += m; tot[2] += h
        L.append(f'| {d+1} | {chunk[0][0]:.1f}–{chunk[-1][0]:.1f} | {len(chunk):,} | {v:,} | {100*m/v:.3f} | {100*h/v:.3f} | {100*(m+h)/v:.3f} |')
    L.append(f'| all | | {n:,} | {tot[0]:,} | {100*tot[1]/tot[0]:.3f} | {100*tot[2]/tot[0]:.3f} | {100*(tot[1]+tot[2])/tot[0]:.3f} |')
open(out, 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
