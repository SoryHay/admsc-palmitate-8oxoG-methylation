#!/usr/bin/env python3
"""5mC and 5hmC within ribosomal DNA interval sets (Table S5).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, bisect, collections, time
P = '/data/8oxo_project/AdMSC_paper'
BEDS = {'rrna': f'{P}/sessions/2026-08-26b_manuscript/rdna_oxo/rdna_T2T.bed', 'acro5': f'{P}/sessions/2026-08-26b_manuscript/rdna_oxo/rdna_acro5.bed'}
PILE = {'NT': f'{P}/sessions/2026-09-03/p1_5hmC/NT_cpg_thr07.bed', 'PA': f'{P}/sessions/2026-09-03/p1_5hmC/PA_cpg_thr07.bed'}
OUT = f'{P}/sessions/2026-09-06/results/p2_rdna_5hmC'
GENOME = {'NT': (66.325, 2.795), 'PA': (64.925, 4.190)}
t0 = time.time()

def load_intervals(path):
    iv = collections.defaultdict(list)
    for line in open(path):
        if not line.strip() or line.startswith('#'): continue
        c = line.split('\t'); iv[c[0]].append((int(c[1]), int(c[2])))
    for k in iv: iv[k].sort()
    starts = {k: [a for a, b in v] for k, v in iv.items()}
    return iv, starts, sum(len(v) for v in iv.values())
def inside(iv, starts, chrom, pos):
    if chrom not in iv: return False
    i = bisect.bisect_right(starts[chrom], pos) - 1
    return i >= 0 and iv[chrom][i][0] <= pos < iv[chrom][i][1]

def summarise(rows):
    tot = collections.Counter(); seen = set()
    for c in rows:
        key = (c[0], c[1], c[5]); valid, nmod, code = int(c[9]), int(c[11]), c[3]
        if key not in seen: seen.add(key); tot['valid'] += valid; tot['positions'] += 1
        if code == 'm': tot['m'] += nmod
        elif code == 'h': tot['h'] += nmod
    return tot

sets = {name: load_intervals(path) for name, path in BEDS.items()}
kept = {name: {arm: [] for arm in PILE} for name in sets}
nrows = {}
for arm, path in PILE.items():
    n = 0
    with open(path) as f:
        for line in f:
            c = line.rstrip('\n').split('\t'); n += 1
            pos = int(c[1])
            for name, (iv, starts, _) in sets.items():
                if inside(iv, starts, c[0], pos): kept[name][arm].append(c)
    nrows[arm] = n
    print(f'{arm}: {n:,} pileup rows read; kept rrna {len(kept["rrna"][arm]):,}, acro5 {len(kept["acro5"][arm]):,}  ({time.time()-t0:.0f} s)', file=sys.stderr)
for name in sets:
    for arm in PILE:
        with open(f'{OUT}/{arm}_{name}.bed', 'w') as g:
            for c in kept[name][arm]: g.write('\t'.join(c) + '\n')

L = ['# P-2 — 5mC and 5hmC in ribosomal DNA, combined model 5mCG_5hmCG@v0, both arms', '',
     f'Run 2026-09-07 by `sessions/2026-09-06/scripts/p2_rdna_5hmc.py`. Inputs: the P-1 whole-genome pileups '
     f'(`sessions/2026-09-03/p1_5hmC/{{NT,PA}}_cpg_thr07.bed`, {nrows["NT"]:,} and {nrows["PA"]:,} rows; the calls behind L20–L27), '
     'filtered to the rDNA intervals of `sessions/2026-08-26b_manuscript/rdna_oxo/` (the same intervals as L14–L19); summation as in '
     '`p1_step8_delta.py`. No basecalling and no new modkit run.', '',
     '## Criterion, fixed before the run (edit sheet Part C, P-2)', '',
     'On Δ5hmC(rDNA) = PA − NT in percentage points of valid CpG calls, against the genome-wide Δ5hmC = **+1.395 pp** (L20–L27):',
     '- **absent**: Δ ≤ 0, or 0 < Δ < 0.70 pp (less than half the genome-wide gain);',
     '- **present**: 0.70 ≤ Δ ≤ 2.10 pp (within ±50 % of the genome-wide gain);',
     '- **larger**: Δ > 2.10 pp.', '',
     'Genome-wide reference (L20–L27): NT 5mC 66.325 %, 5hmC 2.795 %; PA 5mC 64.925 %, 5hmC 4.190 %; Δ5mC −1.400 pp, Δ5hmC +1.395 pp.', '']
verd = {}
for name, label in (('rrna', 'rRNA class, 605 intervals'), ('acro5', 'five acrocentric 45S arrays, 216 segments')):
    L += [f'## {label}', '', '| arm | CpG positions | valid calls | 5mC calls | 5hmC calls | 5mC % | 5hmC % | genome 5mC % | genome 5hmC % |', '|---|---|---|---|---|---|---|---|---|']
    fr = {}
    for arm in ('NT', 'PA'):
        t = summarise(kept[name][arm]); fm, fh = 100 * t['m'] / t['valid'], 100 * t['h'] / t['valid']; fr[arm] = (fm, fh, t)
        L.append(f"| {arm} | {t['positions']:,} | {t['valid']:,} | {t['m']:,} | {t['h']:,} | {fm:.3f} | {fh:.3f} | {GENOME[arm][0]:.3f} | {GENOME[arm][1]:.3f} |")
    dm = fr['PA'][0] - fr['NT'][0]; dh = fr['PA'][1] - fr['NT'][1]
    v = 'absent' if dh <= 0 or dh < 0.70 else ('present' if dh <= 2.10 else 'larger')
    verd[name] = (dm, dh, v, fr)
    L += ['', f'Δ5mC(rDNA) = **{dm:+.3f} pp** (ratio PA/NT {fr["PA"][0]/fr["NT"][0]:.4f}); Δ5hmC(rDNA) = **{dh:+.3f} pp** (ratio {fr["PA"][1]/fr["NT"][1]:.4f}); '
          f'Δ5hmC / |Δ5mC| = {dh/abs(dm) if dm else float("nan"):+.3f}; Δ5hmC(rDNA) / Δ5hmC(genome) = **{dh/1.395:.2f}**.',
          f'Share of the genome-wide valid calls that lie in these intervals: NT {100*fr["NT"][2]["valid"]/5274382:.2f} %, PA {100*fr["PA"][2]["valid"]/5890585:.2f} %.',
          '', f'**Verdict by the criterion: conversion {v.upper()} in ribosomal DNA.**', '']
L += ['## What would overturn it', '',
      'The intervals are the RepeatMasker rRNA class and the acrocentric arrays on T2T; reads mapping into the arrays carry the usual '
      'multi-copy ambiguity, which the pooled fraction tolerates but a per-copy claim would not. The model’s 5hmC state absorbs 9–11 % of '
      'unmodified CpG on the amplified control, so the absolute levels are not claimed; only the arm-to-arm difference within '
      'the model is read, as for L20–L27. A different rDNA definition (the 45S transcription unit only, or the intergenic spacer only) is '
      'the test that has not been run.', '',
      '']
open(f'{OUT}/RESULT_P2_rDNA_5hmC_2026-09-07.md', 'w').write('\n'.join(L) + '\n')
for name, (dm, dh, v, fr) in verd.items():
    print(f'{name}: NT 5mC {fr["NT"][0]:.3f} 5hmC {fr["NT"][1]:.3f} | PA 5mC {fr["PA"][0]:.3f} 5hmC {fr["PA"][1]:.3f} | d5mC {dm:+.3f} d5hmC {dh:+.3f} -> {v}')
print('wrote', f'{OUT}/RESULT_P2_rDNA_5hmC_2026-09-07.md', f'({time.time()-t0:.0f} s)')
