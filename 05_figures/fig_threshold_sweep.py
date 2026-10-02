#!/usr/bin/env python3
"""Figure S9: pooled ratio against the call threshold.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

P = '/data/8oxo_project/AdMSC_DMR'
OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-11/figures'
RES = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-11/results'
os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
plt.rcParams.update({'font.size': 8.5, 'axes.labelsize': 9.5, 'xtick.labelsize': 8.5,
                     'ytick.labelsize': 8.5, 'legend.fontsize': 7.8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white'})
GENOME_WIDE = 0.89845

CLASSES = [('islands', 'CpG island', '#0072B2'),
           ('shores', 'shore', '#009E73'),
           ('promoters_tss1kb', 'promoter', '#D55E00'),
           ('genebodies', 'gene body', '#CC79A7')]
THR = [('0.70', 'layer1_cpg_{}.tsv'), ('0.90', 'layer1_sweep_{}_thr090.tsv'),
       ('0.98', 'layer1_sweep_{}_thr098.tsv')]

def load(path):
    """region key -> (NT_valid, NT_mod, PA_valid, PA_mod)"""
    d = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        try:
            a, b = int(r['NT_valid']), int(r['NT_mod'])
            c, e = int(r['PA_valid']), int(r['PA_mod'])
        except (ValueError, KeyError):
            continue
        if a == 0 or c == 0:
            continue
        d[(r['chrom'], r['start'], r['end'])] = (a, b, c, e)
    return d

def pooled_over(d, keys=None):
    ks = d.keys() if keys is None else [k for k in keys if k in d]
    ntv = ntm = pav = pam = 0; n = 0
    for k in ks:
        a, b, c, e = d[k]
        ntv += a; ntm += b; pav += c; pam += e; n += 1
    return (pam / pav) / (ntm / ntv), ntv + pav, n

data, data_m = {}, {}
out = open(f'{RES}/threshold_sweep.tsv', 'w')
out.write("class\tthreshold\tsource_table\tregion_set\tn_regions\ttotal_valid_calls\t"
          "pooled_ratio\n")
for key, label, col in CLASSES:
    tabs = {thr: tmpl.format(key) for thr, tmpl in THR}
    loaded = {}
    for thr, tab in tabs.items():
        path = f'{P}/results/{tab}'
        if os.path.exists(path):
            loaded[thr] = load(path)
        else:
            print(f"  MISSING {tab}")
    common = set.intersection(*[set(d) for d in loaded.values()]) if loaded else set()
    data[key], data_m[key] = [], []
    for thr, _ in THR:
        if thr not in loaded:
            data[key].append(None); data_m[key].append(None); continue
        ra, oa, na = pooled_over(loaded[thr])
        rm, om, nm = pooled_over(loaded[thr], common)
        data[key].append((float(thr), ra, oa, na))
        data_m[key].append((float(thr), rm, om, nm))
        out.write(f"{label}\t{thr}\t{tabs[thr]}\tall\t{na}\t{oa}\t{ra:.5f}\n")
        out.write(f"{label}\t{thr}\t{tabs[thr]}\tmatched_all_three\t{nm}\t{om}\t{rm:.5f}\n")
        print(f"  {label:<11} thr {thr}   all n={na:>7,} ratio {ra:.4f}   |   "
              f"matched n={nm:>7,} ratio {rm:.4f}   calls {om:>12,}")
out.close()
print(f"\nwrote {RES}/threshold_sweep.tsv")

fig, (ax, axb) = plt.subplots(2, 1, figsize=(6.4, 5.4), sharex=True,
                              gridspec_kw={'height_ratios': [2.0, 1.0], 'hspace': 0.16})
ax.axhline(GENOME_WIDE, color='#444444', lw=1.0, ls='--', zorder=1)
ax.text(0.985, GENOME_WIDE, f' genome-wide {GENOME_WIDE}\n at threshold 0.70',
        ha='right', va='bottom', fontsize=7.0, color='#444444')
for key, label, col in CLASSES:
    pts = [p for p in data[key] if p]
    ptm = [p for p in data_m[key] if p]
    if not pts: continue
    x = [p[0] for p in pts]; y = [p[1] for p in pts]
    ax.plot([p[0] for p in ptm], [p[1] for p in ptm], '-o', color=col, lw=1.6, ms=6,
            mec='white', mew=0.8, label=f'{label} (matched regions)')
    ax.plot(x, y, ':', color=col, lw=1.0, alpha=0.55)
    ax.annotate(f'{ptm[0][1]:.4f}', (ptm[0][0], ptm[0][1]), textcoords='offset points',
                xytext=(-6, 0), ha='right', va='center', fontsize=7.0, color=col)
    ax.annotate(f'{ptm[-1][1]:.4f}', (ptm[-1][0], ptm[-1][1]), textcoords='offset points',
                xytext=(7, 0), ha='left', va='center', fontsize=7.0, color=col)
ax.set_ylabel('pooled PA / NT ratio')
ax.set_xlim(0.63, 1.05)
ax.legend(frameon=False, ncol=2, loc='upper left', fontsize=7.2)
ax.set_title('A   The PA/NT ratio against the caller threshold\n'
             '     solid = regions present at ALL THREE thresholds; dotted = all regions at that '
             'threshold (confounded by dropout)',
             loc='left', fontsize=9)

for key, label, col in CLASSES:
    pts = [p for p in data[key] if p]
    if not pts: continue
    base = pts[0][2]
    axb.plot([p[0] for p in pts], [100 * p[2] / base for p in pts], '-o',
             color=col, lw=1.3, ms=5, mec='white', mew=0.8)
axb.set_ylabel('calls retained (% of 0.70)')
axb.set_xlabel('caller threshold  (P(5mC) required for a call to count)')
axb.set_xticks([0.70, 0.90, 0.98])
axb.set_ylim(0, 108)
axb.axhline(100, color='#cccccc', lw=0.7, ls=':')
axb.set_title('B   How much data each threshold keeps — a stable ratio on 10 % of the calls is a '
              'weaker statement', loc='left', fontsize=9)

fig.text(0.005, -0.02,
         'Pooled ratios only. No q-value anywhere. This is a sensitivity display: it '
         'shows whether the ratio depends on the operating point,\nnot whether any class differs from '
         'any other. n = 1 per arm, condition confounded with flow cell — direction and rank only.  '
         'Values: results/threshold_sweep.tsv',
         fontsize=6.4, color='#444444')
for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT}/Fig_threshold_sweep.{ext}', dpi=600, bbox_inches='tight', facecolor='white')
print(f"-> {OUT}/Fig_threshold_sweep.png")
