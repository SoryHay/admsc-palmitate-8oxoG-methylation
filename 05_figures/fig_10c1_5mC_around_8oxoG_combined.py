#!/usr/bin/env python3
"""Figure 10.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv
import os
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle

HERE = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-26b_manuscript'
OUT8 = '/data/8oxo_project/AdMSC_paper/sessions/2026-10-02'
import os, sys; sys.path.insert(0, '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03/scripts')
from panel_export import export_panels
SRC = f'{HERE}/results/d3_profile_bin50_COPY_from_2026-08-14.tsv'
SURFACE, INK, INK2, MUTED = '#ffffff', '#0b0b0b', '#52514e', '#8a8983'
BLUE, ORANGE, AQUA, RED = '#2a78d6', '#eb6834', '#1baf7a', '#a8462a'
READ, EDGE = '#e6e5e0', '#cdccc6'
BASELINE = {'NT': 0.3297, 'PA': 0.3094}
ANCHORS = {'NT': 1189, 'PA': 1217}

fig = plt.figure(figsize=(12.6, 14.2))
fig.patch.set_facecolor(SURFACE)
outer = fig.add_gridspec(3, 1, height_ratios=[1.18, 0.82, 1.90], hspace=0.26,
                         top=0.975, bottom=0.045, left=0.062, right=0.978)

def schem(i, title, sub):
    ax = fig.add_subplot(outer[i]); ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
    ax.text(0, 100, title.split(')')[0] + ')', fontsize=13, color=INK, va='top', fontweight='semibold')
    return ax

def band(ax, x0, x1, y, h=5.0):
    ax.add_patch(FancyBboxPatch((x0, y - h/2), x1 - x0, h,
                                boxstyle='round,pad=0,rounding_size=1.1',
                                facecolor=READ, edgecolor=EDGE, lw=1))

ax = schem(0, '(a)   What is read off one molecule',
           'Nanopore reads the whole molecule once; two independent callers are then run over that same '
           'read. The values shown are illustrative.')
X0, X1, Y = 24, 97, 50
band(ax, X0, X1, Y)
ax.text((X0 + X1) / 2, Y, 'one read · a single DNA molecule, ~1,800 bp',
        ha='center', va='center', fontsize=9, color=MUTED, style='italic')
ax.text(X0, 83.5, 'Guppy + Rerio  →  every CpG on this read gets  P(5mC), a probability in [0, 1]',
        fontsize=9.6, color='#12805a')
for x, p in [(30, 0.88), (39, 0.14), (46, 0.62), (55, 0.09), (64, 0.71), (72, 0.22), (82, 0.55), (92, 0.41)]:
    ax.plot([x, x], [Y + 2.9, Y + 7.5], color=AQUA, lw=2.3, solid_capstyle='round')
    ax.text(x, 58.3, 'C', ha='center', va='bottom', fontsize=10, color='#12805a', fontweight='bold')
    ax.text(x, 63.0, f'{p:.2f}', ha='center', va='bottom', fontsize=9, color='#12805a')
for xc, lab, yy in ((46, 'd = 210 bp', 71.5), (82, 'd = 780 bp', 78.5)):
    ax.add_patch(FancyArrowPatch((30, yy), (xc, yy), arrowstyle='<->', color=INK2,
                                 lw=1.1, mutation_scale=9))
    ax.text((30 + xc) / 2, yy + 1.3, lab, ha='center', fontsize=8.8, color=INK2)
ax.text(X0, 21, 'esox  →  every G whose 5-mer is on the 110-mer allowlist gets a score in [0, 1]',
        fontsize=9.6, color=INK2)
for x, sc, is_oxo in ((30, 0.87, True), (44, 0.00, False), (58, 0.02, False),
                      (75, 0.00, False), (90, 0.31, False)):
    col = BLUE if is_oxo else MUTED
    ax.plot([x, x], [Y - 2.9, Y - 7.5], color=col, lw=2.3, solid_capstyle='round')
    ax.text(x, 38.5, 'G', ha='center', va='top', fontsize=10, color=col, fontweight='bold')
    ax.text(x, 30.5, f'{sc:.2f}', ha='center', va='top', fontsize=9.5, color=col,
            fontweight='bold' if is_oxo else 'normal')
    if is_oxo:
        ax.add_patch(Circle((x, 40.8), 3.0, fill=False, edgecolor=BLUE, lw=1.9, zorder=5))
ax.text(31, 15, 'score ≥ 0.70  →  8-oxo-dG ANCHOR', ha='center', va='top',
        fontsize=9.4, color=BLUE, fontweight='medium')
ax.text(31, 10, 'feeds the blue curve', ha='center', va='top', fontsize=9.0, color=BLUE)
ax.text(72, 15, 'score < 0.70  →  ORDINARY G ANCHOR', ha='center', va='top',
        fontsize=9.4, color=ORANGE, fontweight='medium')
ax.text(72, 10, 'feeds the orange and green curves', ha='center', va='top', fontsize=9.0, color=ORANGE)
ax.add_patch(Rectangle((28.0, Y - 2.5), 4.0, 5.0, facecolor='#f6d9cd', edgecolor='none', zorder=1))
ax.text(0, 56, '±10 bp\nmasked', fontsize=9.3, color=RED, va='center', fontweight='medium',
        linespacing=1.4)
ax.text(0, 47, 'the 5mC model has no\nrepresentation of 8-oxo-dG,\nso a lesion right beside a\n'
               'cytosine distorts the signal\nit reads', fontsize=8.4, color=RED, va='top',
        linespacing=1.5)
ax.add_patch(FancyArrowPatch((15.5, 55), (27.6, 51.5), arrowstyle='->', color=RED, lw=1.1,
                             mutation_scale=10, connectionstyle='arc3,rad=-0.15'))
ax.text((X0 + X1) / 2, 1.5,
        'Every CpG on this read lying 11–1000 bp from an anchor forms one “anchor–CpG pair”. '
        'Not drawn to scale.', ha='center', fontsize=9.5, color=INK)

ax_a = ax

ax = schem(1, '(b)   Which anchors feed which curve',
           'All three curves plot the same quantity — mean P(5mC). They differ only in which anchors '
           'they are measured around.')
for yy, name in ((64, 'a molecule that\nCARRIES a lesion'), (26, 'a molecule with\nNO lesion at all')):
    band(ax, 20, 58, yy)
    ax.text(18.5, yy, name, ha='right', va='center', fontsize=9.3, color=INK2, linespacing=1.5)
for x, oxo in ((28, True), (38, False), (48, False)):
    col = BLUE if oxo else MUTED
    ax.plot([x, x], [64 - 2.9, 64 - 6.6], color=col, lw=2.3, solid_capstyle='round')
    if oxo:
        ax.add_patch(Circle((x, 60.0), 2.7, fill=False, edgecolor=BLUE, lw=1.9, zorder=5))
for x in (28, 38, 48):
    ax.plot([x, x], [26 - 2.9, 26 - 6.6], color=MUTED, lw=2.3, solid_capstyle='round')
for (xa, ya), (xt, yt), col, lab in (
        ((28, 56.8), (62, 78), BLUE,   'BLUE — measured around the oxidised G'),
        ((43, 56.8), (62, 60), ORANGE, 'ORANGE — measured around the ordinary Gs\nof that same molecule'),
        ((38, 18.8), (62, 24), AQUA,   'GREEN — measured around the ordinary Gs of\nmolecules carrying no lesion at all')):
    ax.add_patch(FancyArrowPatch((xa, ya), (xt - 1.5, yt), arrowstyle='->', color=col, lw=1.6,
                                 mutation_scale=12, connectionstyle='arc3,rad=-0.18'))
    ax.text(xt, yt, lab, va='center', fontsize=9.3, color=col, linespacing=1.5)
ax.text(50, 2,
        'The ordinary Gs of a lesion-carrying molecule feed ORANGE and never GREEN — which is why '
        'ORANGE against GREEN measures what KIND\nof molecule carries lesions, with no reference to '
        'distance at all.', ha='center', fontsize=9.3, color=INK, linespacing=1.6, va='top')

ax_b = ax

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
        return ([(r['dist_lo'] + r['dist_hi']) / 2 for r in c], [r['mean'] for r in c], [int(r['n']) for r in c])
    s = sorted((r for r in rows if r['arm'] == arm and r['class'] == cls),
               key=lambda r: int(r['dist_lo']))
    return ([(int(r['dist_lo']) + int(r['dist_hi'])) / 2 for r in s],
            [float(r['mean_prediction_score']) for r in s],
            [int(r['n_pairs']) for r in s])

inner = outer[2].subgridspec(2, 2, height_ratios=[3.3, 1.0], hspace=0.13, wspace=0.13)
SERIES = [('oxo', 'oxidised G', BLUE, 'o'),
          ('canon_same', 'ordinary G,\nsame molecules', ORANGE, 's'),
          ('canon_free', 'ordinary G,\nundamaged molecules', AQUA, '^')]
_c_axes = []
for col, arm in enumerate(('NT', 'PA')):
    ax = fig.add_subplot(inner[0, col]); axn = fig.add_subplot(inner[1, col], sharex=ax); _c_axes += [ax, axn]
    for a in (ax, axn):
        a.set_facecolor(SURFACE)
        for s_ in ('top', 'right'):
            a.spines[s_].set_visible(False)
        for s_ in ('left', 'bottom'):
            a.spines[s_].set_color('#d8d7d2')
        a.tick_params(colors=INK2, labelsize=8.5, length=3)
        a.grid(axis='y', color='#ededea', lw=0.8); a.set_axisbelow(True)
    ax.axhline(BASELINE[arm], color='#9c9b95', lw=1.2, ls=(0, (5, 4)), zorder=1)
    ax.text(18, BASELINE[arm] + 0.0035, 'all-CpG mean', ha='left', va='bottom',
            fontsize=8.5, color=INK2)
    for cls, lab, colour, marker in SERIES:
        x, y, _ = ser(arm, cls)
        ax.plot(x, y, color=colour, lw=2, marker=marker, ms=4.5, mew=0, label=lab, zorder=3)
        if col == 0:
            ax.annotate(lab, xy=(x[-1], y[-1]), xytext=(8, 0), textcoords='offset points',
                        va='center', fontsize=9, color=colour, fontweight='medium')
    ax.set_ylim(0.20, 0.40); ax.set_xlim(0, 1265)
    ax.text(0.012, 0.975, f'{arm} · {ANCHORS[arm]:,} 8-oxo-dG anchors', transform=ax.transAxes,
            fontsize=10.5, color=INK, va='top')
    ax.tick_params(labelbottom=False)
    if col == 0:
        _first = ax
    xo, _, n = ser(arm, 'oxo')
    axn.bar(xo, n, width=42, color='#c9c8c2', zorder=2)
    axn.set_ylim(0, 2600); axn.set_xlim(0, 1265)
    axn.set_xticks([0, 250, 500, 750, 1000])
    axn.set_xlabel('distance from the anchor on the reference (bp)', fontsize=9.5, color=INK)
    axn.annotate(f'{n[0]:,} → {n[-1]:,} pairs', xy=(1250, 2400), ha='right', fontsize=8.5,
                 color=INK2, va='top')
    if col == 0:
        ax.set_ylabel('mean 5mC prediction score\nof the CpGs at that distance', fontsize=9.5, color=INK)
        axn.set_ylabel('anchor–CpG pairs\nper 50 bp bin', fontsize=9, color=INK)
        ax.legend(loc='lower center', frameon=False, fontsize=8.5, labelcolor=INK2,
                  handlelength=1.6, borderpad=0.2)
    else:
        ax.tick_params(labelleft=False); axn.tick_params(labelleft=False)

_p = _first.get_position()
fig.text(0.062, _p.y1 + 0.040, '(c)', fontsize=13, color=INK, fontweight='semibold', va='top')

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT8}/figures/Fig10-1_5mC_around_8oxoG_combined.{ext}', dpi=DPI, facecolor=SURFACE)
print('wrote figures/Fig_5mC_around_8oxoG_combined.png and .pdf')

