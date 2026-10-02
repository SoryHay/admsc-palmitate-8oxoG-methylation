#!/usr/bin/env python3
"""Figure S6a-d: gene panels, gene by gene.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import sys, csv, os, subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({'font.size': 7.5, 'axes.labelsize': 8.5, 'xtick.labelsize': 7.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white'})
P = '/data/8oxo_project/AdMSC_DMR'
OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-11/figures'
os.makedirs(OUT, exist_ok=True)
R_GLOBAL = 0.8971
MIN_NT = 10.0

PANELS = {
    'oxidative':      ('Oxidative damage, repair and redox', 'oxidative_genes.tsv'),
    'innate_exons':   ('Innate immunity',                    'innate_immunity_genes.tsv'),
    'mito_metabolic': ('Mitochondrial and metabolic',        'mito_metabolic_genes.tsv'),
    'adipokine':      ('Adipokines measured by the assay',   'adipokine_assay_genes.tsv'),
    'senescence':     ('Ageing and early senescence',         'senescence_genes.tsv'),
}

def overlaps(panel):
    """fraction of each gene's span covered by a different gene"""
    bed = f'{P}/features/{panel}_bygene.bed'
    if not os.path.exists(bed):
        bed = f'{P}/features/{panel.replace("_exons", "")}_bygene.bed'
    res = subprocess.run(
        f"bedtools intersect -a <(cut -f1-4 {bed}|sort -k1,1 -k2,2n) "
        f"-b <(sort -k1,1 -k2,2n {P}/features/genebodies.bed|cut -f1-4) -wao",
        shell=True, capture_output=True, text=True, executable='/bin/bash')
    if res.returncode != 0 or not res.stdout.strip():
        raise SystemExit(f"FATAL: bedtools produced nothing for {panel}: {res.stderr[:300]}\n"
                         f"the overlap annotation is required; refusing to draw without it.")
    ov, ln = {}, {}
    for line in res.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 9 or f[3] == f[7]:
            continue
        try:
            n = int(f[8])
        except ValueError:
            continue
        if n > 0:
            ov[f[3]] = ov.get(f[3], 0) + n
            ln[f[3]] = int(f[2]) - int(f[1])
    return {g: ov[g] / ln[g] for g in ov if ln.get(g)}

for panel, (title, famfile) in PANELS.items():
    if len(sys.argv) > 1 and panel not in sys.argv[1:]:
        continue
    tsv = f'{P}/results/layer1_cpg_{panel}_bygene.tsv'
    if not os.path.exists(tsv):
        print('missing', tsv); continue
    fam = {}
    fp = f'{P}/features/{famfile}'
    if os.path.exists(fp):
        for l in open(fp):
            if l.startswith('#'): continue
            c = l.rstrip('\n').split('\t')
            if len(c) >= 2: fam[c[0]] = c[1]
    ovf = overlaps(panel)

    rows = []
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        try:
            nt, pa = float(r['NT_rate']), float(r['PA_rate'])
            if nt <= 0: continue
            rows.append(dict(g=r['name'], fam=fam.get(r['name'], 'other'),
                             nt=100*nt, pa=100*pa, ratio=pa/nt,
                             n=min(int(r['NT_valid']), int(r['PA_valid'])),
                             ov=ovf.get(r['name'], 0.0)))
        except (ValueError, KeyError):
            continue
    if not rows: continue
    fams = sorted({r['fam'] for r in rows})
    grouped = 1 < len(fams) < len(rows) * 0.7
    if grouped:
        ordered = []
        for f in fams:
            ordered += sorted([r for r in rows if r['fam'] == f], key=lambda d: d['ratio'])
    else:
        ordered = sorted(rows, key=lambda d: d['ratio'])
    n = len(ordered)
    weak = np.array([r['nt'] < MIN_NT for r in ordered])
    h = max(2.6, 0.20 * n + 1.2) if n < 30 else max(2.6, 0.155 * n + 1.1)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.9, h),
                                   gridspec_kw={'width_ratios': [1.35, 1]})
    left = 0.30 if not grouped else 0.24
    fig.subplots_adjust(wspace=0.5, left=left, right=0.94, top=1-0.62/h, bottom=0.58/h)

    y = np.arange(n)[::-1]

    vals = np.array([r['ratio'] for r in ordered])
    norm = TwoSlopeNorm(vmin=min(0.3, vals.min()), vcenter=R_GLOBAL, vmax=max(1.8, vals.max()))
    cols = [(0.85, 0.85, 0.85, 1.0) if w else plt.cm.RdBu(norm(v))
            for v, w in zip(vals, weak)]
    axL.barh(y, vals - R_GLOBAL, left=R_GLOBAL, height=0.72, color=cols, edgecolor='none')
    axL.axvline(R_GLOBAL, color='k', lw=0.9)
    axL.text(R_GLOBAL, n + 0.4, f'gene-body\nbackground {R_GLOBAL:.4f}',
             ha='center', va='bottom', fontsize=6)
    labels = []
    for r in ordered:
        mark = '●' if r['ov'] >= 0.5 else ('○' if r['ov'] > 0 else ' ')
        dag = '†' if r['nt'] < MIN_NT else ' '
        extra = '' if (grouped or r['fam'] in ('other', r['g'])) else f"  ({r['fam']})"
        labels.append(f"{mark} {r['g']}{dag}{extra}")
    axL.set_yticks(y); axL.set_yticklabels(labels)
    for tick, w in zip(axL.get_yticklabels(), weak):
        if w: tick.set_color('#a0a0a0')
    axL.set_ylim(-0.8, n + 1.6)
    axL.set_xlabel('PA / NT ratio')
    axL.set_title(f'{title}\nleft of line = lost more than background', loc='left', fontsize=8)

    if grouped:
        i = 0
        for f in fams:
            k = sum(1 for r in ordered if r['fam'] == f)
            axL.text(-0.30, n - (i + k/2) - 0.5, f.replace('_', ' '),
                     transform=axL.get_yaxis_transform(), ha='right', va='center',
                     fontsize=6, rotation=90, color='#444444')
            i += k
            if f != fams[-1]:
                axL.axhline(n - i - 0.5, color='#cccccc', lw=0.5)

    axR.barh(y + 0.19, [r['nt'] for r in ordered], height=0.36, color='#0072B2', label='NT')
    axR.barh(y - 0.19, [r['pa'] for r in ordered], height=0.36, color='#D55E00', label='PA')
    axR.axvline(MIN_NT, color='#a0a0a0', ls=':', lw=0.8)
    axR.set_yticks([]); axR.set_ylim(-0.8, n + 1.6)
    axR.set_xlabel('5mCpG (%)')
    axR.legend(frameon=False, ncol=2, loc='upper right', fontsize=6.5)
    axR.set_title('absolute level   (dotted = 10 % NT)', loc='left', fontsize=8)

    fig.text(0.005, 0.004,
             '† NT < 10 %, ratio not interpretable (greyed)   '
             '○ span overlaps another gene   ● overlap >= 50 % of span\n'
             'NO significance marker: no per-gene test is made. '
             'Effect size only - no gene here is claimed to differ from background.',
             fontsize=6, color='#444444')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{OUT}/Fig_panel_{panel}_noq.{ext}', dpi=600, bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)
    print(f'-> Fig_panel_{panel}_noq  ({n} genes, {int(weak.sum())} greyed NT<10%, '
          f'{sum(1 for r in ordered if r["ov"]>0)} overlapping)')
