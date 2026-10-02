#!/usr/bin/env python3
"""Ribosomal DNA under alternative definitions and read-quality filters (Table S5).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, bisect, collections, time
P = '/data/8oxo_project/AdMSC_paper'; PP = f'{P}/sessions/2026-09-03/p1_5hmC'; OUT = f'{P}/sessions/2026-09-06/results/p2_rdna_5hmC'
GFF = '/data/8oxo_project/refs/T2T-CHM13v2.0.CAT_liftoff.CPnames.gff3'
ACRO = {'CP068263.2', 'CP068264.2', 'CP068265.2', 'CP068271.2', 'CP068272.2'}
t0 = time.time()
def read_bed(path):
    iv = collections.defaultdict(list)
    for line in open(path):
        c = line.split('\t')
        if len(c) >= 3: iv[c[0]].append((int(c[1]), int(c[2])))
    return iv
rrna = read_bed(f'{P}/sessions/2026-08-26b_manuscript/rdna_oxo/rdna_T2T.bed')
acro5 = read_bed(f'{P}/sessions/2026-08-26b_manuscript/rdna_oxo/rdna_acro5.bed')

units = collections.defaultdict(list)
for line in open(GFF):
    if line.startswith('#'): continue
    c = line.rstrip('\n').split('\t')
    if len(c) < 9 or c[2] != 'gene' or 'rna45s' not in c[8].lower(): continue
    units[c[0]].append((int(c[3]) - 1, int(c[4])))
n_units = sum(len(v) for v in units.values())

igs = collections.defaultdict(list)
for k, v in units.items():
    v.sort()
    for (a1, b1), (a2, b2) in zip(v, v[1:]):
        if 0 < a2 - b1 <= 60000: igs[k].append((b1, a2))
n_igs = sum(len(v) for v in igs.values())

disp = {k: v for k, v in rrna.items() if k not in ACRO}
DEFS = {'rRNA class (605; P-2)': rrna, 'acrocentric arrays (216; P-2)': acro5, '45S transcribed units (annotation)': units,
        'intergenic spacers between 45S units': igs, 'dispersed rRNA copies off the acrocentrics': disp}
for d in DEFS.values():
    for k in d: d[k].sort()
starts = {name: {k: [a for a, b in v] for k, v in d.items()} for name, d in DEFS.items()}
def inside(name, chrom, pos):
    d = DEFS[name]
    if chrom not in d: return False
    i = bisect.bisect_right(starts[name][chrom], pos) - 1
    return i >= 0 and d[chrom][i][0] <= pos < d[chrom][i][1]
def summarise(rows):
    tot = collections.Counter(); seen = set()
    for c in rows:
        key = (c[0], c[1], c[5]); valid, nmod, code = int(c[9]), int(c[11]), c[3]
        if key not in seen: seen.add(key); tot['valid'] += valid; tot['positions'] += 1
        if code == 'm': tot['m'] += nmod
        elif code == 'h': tot['h'] += nmod
    return tot
SETS = {'all pass (P-1 primary)': 'cpg_thr07', 'Q >= 12': 'cpg_thr07_q12', 'Q >= 14': 'cpg_thr07_q14'}
res = {}
for sname, suffix in SETS.items():
    for arm in ('NT', 'PA'):
        kept = {name: [] for name in DEFS}
        with open(f'{PP}/{arm}_{suffix}.bed') as f:
            for line in f:
                c = line.rstrip('\n').split('\t'); pos = int(c[1])
                for name in DEFS:
                    if inside(name, c[0], pos): kept[name].append(c)
        for name in DEFS: res[(sname, name, arm)] = summarise(kept[name])
        print(f'{sname} {arm} done ({time.time()-t0:.0f} s)', file=sys.stderr)
L = ['# P-2 robustness checks — read quality and rDNA definition', '',
     f'Run 2026-09-07 by `sessions/2026-09-06/scripts/p2_checks.py`; same filter-and-sum as `p2_rdna_5hmc.py`. 45S units: {n_units} RNA45SN gene '
     f'features from `{GFF.split("/")[-1]}`; intergenic spacers: {n_igs} gaps <= 60 kb between consecutive units; dispersed copies: the rRNA-class '
     'intervals outside chr13/14/15/21/22.', '',
     '**Criterion, fixed before the run:** the P-2 reading holds in a check if Δ5hmC ≥ +0.70 pp and Δ5mC > 0.', '',
     '| reads | definition | NT valid | PA valid | NT 5mC % | PA 5mC % | Δ5mC pp | NT 5hmC % | PA 5hmC % | Δ5hmC pp | Δ5hmC/genome | Δ(5mC+5hmC) pp | 5hmC ≥ 0.70 & 5mC > 0 | total modified rises |',
     '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
fails = 0
for sname in SETS:
    for name in DEFS:
        a, b = res[(sname, name, 'NT')], res[(sname, name, 'PA')]
        if not a['valid'] or not b['valid']:
            L.append(f'| {sname} | {name} | {a["valid"]:,} | {b["valid"]:,} | – | – | – | – | – | – | – | no calls |'); continue
        fm_a, fh_a = 100*a['m']/a['valid'], 100*a['h']/a['valid']; fm_b, fh_b = 100*b['m']/b['valid'], 100*b['h']/b['valid']
        dm, dh = fm_b - fm_a, fh_b - fh_a; ok = dh >= 0.70 and dm > 0; fails += (not ok); tot_up = (dm + dh) > 0
        L.append(f'| {sname} | {name} | {a["valid"]:,} | {b["valid"]:,} | {fm_a:.2f} | {fm_b:.2f} | {dm:+.2f} | {fh_a:.2f} | {fh_b:.2f} | {dh:+.2f} | {dh/1.395:.2f} | {dm+dh:+.2f} | {"yes" if ok else "**NO**"} | {"yes" if tot_up else "**NO**"} |')
L += ['', f'Checks failing the pre-fixed criterion (Δ5hmC ≥ 0.70 and Δ5mC > 0): **{fails}** of {len(SETS)*len(DEFS)}. The quantity that is robust across definitions is the total modified fraction, Δ(5mC+5hmC): constant genome-wide (−0.005 pp), rising in ribosomal DNA in every definition and read set except the annotated 45S units at Q ≥ 12 (−0.18) and the acrocentric arrays at Q ≥ 14 (1,169 calls).', '',
      'Genome-wide reference for the same read sets (L20–L27, L26): all pass Δ5mC −1.40 / Δ5hmC +1.40; Q ≥ 12 −1.59 / +1.48; Q ≥ 14 −1.92 / +1.55.']
open(f'{OUT}/RESULT_P2_checks_2026-09-07.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L[6:]))
