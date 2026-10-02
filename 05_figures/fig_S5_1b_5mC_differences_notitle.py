#!/usr/bin/env python3
"""Figure S5.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-26b_manuscript'
SRC = f'{HERE}/results/d3_profile_bin50_COPY_from_2026-08-14.tsv'
ANCHORS = {'NT': 1189, 'PA': 1217}
SURFACE, INK, INK2 = '#fcfcfb', '#0b0b0b', '#52514e'
C_PROX, C_COMP = '#2a78d6', '#eb6834'

rows = [r for r in csv.DictReader(open(SRC), delimiter='\t') if r['cut'] == '0.7']

def _canon_same(arm):
    """canon_same per bin = (mean_all*n_all - mean_free*n_free) / (n_all - n_free)  (exact; d4 derivation)."""
    def _rows(cls):
        return {int(r['dist_lo']): r for r in rows if r['arm'] == arm and r['class'] == cls}
    a, f = _rows('canon_all'), _rows('canon_free')
    assert a.keys() == f.keys()
    out = []
    for lo in sorted(a):
        na, nf = float(a[lo]['n_pairs']), float(f[lo]['n_pairs'])
        s = float(a[lo]['mean_prediction_score']) * na - float(f[lo]['mean_prediction_score']) * nf
        out.append(dict(arm=arm, cls='canon_same', dist_lo=lo, dist_hi=int(a[lo]['dist_hi']), n=na - nf, mean=s / (na - nf)))
    return out
def ser(arm, cls):
    if cls == 'canon_same':
        c = _canon_same(arm)
        return ([(r['dist_lo'] + r['dist_hi']) / 2 for r in c], [r['mean'] for r in c])
    s = sorted((r for r in rows if r['arm'] == arm and r['class'] == cls),
               key=lambda r: int(r['dist_lo']))
    return ([(int(r['dist_lo']) + int(r['dist_hi'])) / 2 for r in s],
            [float(r['mean_prediction_score']) for r in s])

fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
fig.patch.set_facecolor(SURFACE)

for col, arm in enumerate(('NT', 'PA')):
    ax = axes[col]
    ax.set_facecolor(SURFACE)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=9, length=3)
    ax.grid(axis='y', color='#ededea', lw=0.8)
    ax.set_axisbelow(True)

    x, oxo = ser(arm, 'oxo')
    _, same = ser(arm, 'canon_same')
    _, free = ser(arm, 'canon_free')
    prox = [a - b for a, b in zip(oxo, same)]
    comp = [b - c for b, c in zip(same, free)]

    ax.axhline(0, color='#6f6e69', lw=1.4, zorder=2)
    ax.plot(x, prox, color=C_PROX, lw=2, marker='o', ms=4.5, mew=0, zorder=3,
            label='proximity:  oxidised G  −  ordinary G, same molecules')
    ax.plot(x, comp, color=C_COMP, lw=2, marker='s', ms=4.5, mew=0, zorder=3,
            label='composition:  ordinary G, same molecules  −  undamaged molecules')

    ax.set_ylim(-0.095, 0.02)
    ax.set_xlim(0, 1150)
    ax.set_xticks([0, 250, 500, 750, 1000])
    ax.set_title(f'{arm} · {ANCHORS[arm]:,} 8-oxo-dG anchors',
                 fontsize=11.5, color=INK, loc='left', pad=8)
    ax.set_xlabel('distance from the anchor on the reference (bp)', fontsize=10, color=INK)
    if col == 0:
        ax.set_ylabel('difference in mean 5mC prediction score', fontsize=10, color=INK)
        ax.annotate('no difference', xy=(1140, 0.004), ha='right', fontsize=8.5, color=INK2)
        ax.annotate('less methylation\nin the first term', xy=(1140, -0.086), ha='right',
                    fontsize=8.5, color=INK2, va='center')
    if col == 0:
        for y, c, lab in ((prox[-1], C_PROX, 'proximity'), (comp[-1], C_COMP, 'composition')):
            ax.annotate(lab, xy=(x[-1], y), xytext=(8, 0), textcoords='offset points',
                        va='center', fontsize=9, color=c, fontweight='medium')

h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=1, frameon=False, fontsize=9.5,
           labelcolor=INK2, bbox_to_anchor=(0.5, -0.015), handlelength=1.8)

fig.text(0.055, 0.925,
         'each curve IS one of the two comparisons · flat cut 0.70 · ±10 bp crosstalk mask · '
         'one library per arm',
         fontsize=10, color=INK2, ha='left')
fig.subplots_adjust(top=0.83, bottom=0.255, left=0.085, right=0.965, wspace=0.10)

for ext in ('png', 'pdf'):
    fig.savefig(f'/data/8oxo_project/AdMSC_paper/sessions/2026-10-02/figures/FigS5-1b_5mC_around_8oxoG_differences_notitle.{ext}', dpi=200, facecolor=SURFACE, bbox_inches='tight', pad_inches=0.15)
print('wrote figures/Fig_5mC_around_8oxoG_differences.png and .pdf')
