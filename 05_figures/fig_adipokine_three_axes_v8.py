#!/usr/bin/env python3
"""Figure 12.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os, math
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

S = '/data/8oxo_project/AdMSC_paper/sessions'
HERE = f'{S}/2026-08-26b_manuscript'
OUT8 = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03'
import os, sys; os.makedirs(f'{OUT8}/figures/for_docx', exist_ok=True); sys.path.insert(0, f'{OUT8}/scripts')
from panel_export import export_panels
SURFACE, INK, INK2, MUTED = '#ffffff', '#0b0b0b', '#52514e', '#8a8983'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
PROM_BG = 0.9130

fig = plt.figure(figsize=(12.6, 15.6)); fig.patch.set_facecolor(SURFACE)

gs = fig.add_gridspec(3, 2, height_ratios=[3.97, 3.86, 4.45], width_ratios=[1.0, 1.0],
                      left=0.04, right=0.975, top=0.965, bottom=0.03, hspace=0.28, wspace=0.18)

def letter(ax, s, x=0.0, y=1.03, title='', gap=0.045):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=13, color=INK, fontweight='semibold', va='bottom')
    ax.text(x + gap, y, title, transform=ax.transAxes, fontsize=10, color=INK2, va='bottom')

def embed(ax, path, crop=None):
    ax.set_facecolor(SURFACE); ax.axis('off')
    img = mpimg.imread(path)
    if crop:
        h, w = img.shape[:2]; y0, y1, x0, x1 = crop
        img = img[int(h*y0):int(h*y1), int(w*x0):int(w*x1)]
    ax.imshow(img, interpolation='lanczos')

axa = fig.add_subplot(gs[0, :]); embed(axa, f'{S}/2026-08-20/figures/Fig_adipokine_matrix.png', crop=(0.0, 1.0, 0.0, 1.0))
letter(axa, '(a)')

axb = fig.add_subplot(gs[1, 0]); embed(axb, f'{S}/2026-08-20/figures/Fig_adipokine_protein.png', crop=(0.072, 0.325, 0.555, 1.0))
letter(axb, '(b)', x=-0.02, y=1.04)

GREY = '#9c9b95'
G = [
 ('IL6',    36.21, 28.21,  58, 14.3, 135.2,  663.4,   900.4,  1561.7, 'quant'),
 ('IL1B',   None,  None, None,  2.0,   2.8,   10.3,     9.3,    11.2, 'straddle'),
 ('TNF',    17.59, 12.35, 108, None,  None,    1.3,     1.5,     1.0, 'straddle'),
 ('RETN',   25.64, 18.07, 117, None,  None,    1.0,     0.2,     0.5, 'below'),
 ('CXCL10', 21.88, 17.14,  32, None,  None,    5.4,     7.1,    12.8, 'straddle'),
 ('IL10',   35.53, 37.29,  76, None,  None,    0.2,     0.0,     0.4, 'below'),
 ('RBP4',    1.51,  0.94, 463, None,  None, 18932.3, 22894.3, 19491.8, 'quant'),
 ('CCL2',    2.56,  0.00,  39, None,  None, 1063.8,  1245.8,  1387.1, 'quant'),
]
inner = gs[1, 1].subgridspec(1, 3, wspace=0.10, width_ratios=[1.0, 0.9, 1.15])
axes_c = []
for j in range(3):
    ax = fig.add_subplot(inner[j]); axes_c.append(ax); ax.set_facecolor(SURFACE)
    for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
    for s_ in ('left', 'bottom'): ax.spines[s_].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    ax.grid(axis='x', color='#ededea', lw=0.8); ax.set_axisbelow(True)
    ax.set_ylim(-0.6, len(G) - 0.4); ax.invert_yaxis()
for y, (g, pn, pp, pv, tn, tp, b, p250, p500, st) in enumerate(G):
    ax = axes_c[0]
    if pn is None:
        ax.text(20, y, 'no interval', va='center', ha='center', fontsize=7.4, color=MUTED)
    else:
        rd = pn < 10.0
        ax.plot([pn, pp], [y, y], color='#cdccc6', lw=2, zorder=1)
        ax.plot(pn, y, 'o', ms=7, zorder=3, color=SURFACE if rd else BLUE, markeredgecolor=BLUE, markeredgewidth=1.5)
        ax.plot(pp, y, 'o', ms=7, zorder=3, color=SURFACE if rd else ORANGE, markeredgecolor=ORANGE, markeredgewidth=1.5)
    ax = axes_c[1]
    if tn is None:
        ax.text(20, y, '—', va='center', ha='center', fontsize=9, color=MUTED)
    else:
        ax.plot([tn, tp], [y, y], color='#cdccc6', lw=2, zorder=1)
        ax.plot(tn, y, 's', ms=7, color=BLUE, zorder=3); ax.plot(tp, y, 's', ms=7, color=ORANGE, zorder=3)
    ax = axes_c[2]
    vals = [max(v, 0.15) for v in (b, p250, p500)]
    if st == 'quant':
        ax.plot(vals[0], y, 'o', ms=7, color=BLUE, zorder=3)
        ax.plot(vals[1], y, 'o', ms=7, color=ORANGE, zorder=3)
        ax.plot(vals[2], y, 's', ms=6.5, color=ORANGE, zorder=3)
    else:
        col = AQUA if st == 'straddle' else GREY
        for v, mk in zip(vals, ('^', '^', '^')):
            ax.plot(v, y, mk, ms=8, color=SURFACE, markeredgecolor=col, markeredgewidth=1.6, zorder=3)

for y, (g, pn, pp, pv, tn, tp, b, p250, p500, st) in enumerate(G):
    if pn is not None and pn >= 10.0:
        axes_c[0].text(47.5, y, f'×{pp/pn:.2f}', va='center', ha='right', fontsize=7.2, color=INK2)
    if tn is not None:
        axes_c[1].text(560, y, f'×{tp/tn:.1f}', va='center', ha='right', fontsize=7.2, color=INK2)
    if st == 'quant':
        axes_c[2].text(max(b, p250, p500) * 1.7, y, f'×{p250/b:.2f} · ×{p500/b:.2f}', va='center', ha='left', fontsize=7.0, color=INK2)
axes_c[0].set_yticks(range(len(G))); axes_c[0].set_yticklabels([g[0] for g in G], fontsize=8.5)
for t, g in zip(axes_c[0].get_yticklabels(), G):
    if g[0] == 'IL6': t.set_fontweight('bold')
axes_c[0].set_xlim(0, 48); axes_c[0].set_xlabel('promoter CpG methylation (%)', fontsize=8.6, color=INK)
axes_c[1].tick_params(labelleft=False); axes_c[1].set_xscale('log'); axes_c[1].set_xlim(1, 600)
axes_c[1].set_xlabel('transcript, dPCR (copies/µL)', fontsize=8.6, color=INK)
axes_c[2].tick_params(labelleft=False); axes_c[2].set_xscale('log'); axes_c[2].set_xlim(0.1, 600000)
axes_c[2].set_xlabel('protein at 24 h (pg/mL)', fontsize=8.6, color=INK)
from matplotlib.lines import Line2D
axes_c[2].legend(handles=[
    Line2D([], [], marker='o', color=BLUE, ls='', ms=6, label='untreated / BSA'),
    Line2D([], [], marker='o', color=ORANGE, ls='', ms=6, label='250 µM palmitate'),
    Line2D([], [], marker='s', color=ORANGE, ls='', ms=6, label='500 µM palmitate'),
    Line2D([], [], marker='^', color=SURFACE, markeredgecolor=AQUA, markeredgewidth=1.5, ls='', ms=7, label='straddles its floor'),
    Line2D([], [], marker='^', color=SURFACE, markeredgecolor=GREY, markeredgewidth=1.5, ls='', ms=7, label='at or below its floor'),
    Line2D([], [], marker='o', color=SURFACE, markeredgecolor=BLUE, markeredgewidth=1.5, ls='', ms=6, label='promoter < 10 % — not read'),
], loc='lower left', bbox_to_anchor=(-2.35, -0.50), ncol=3, frameon=False, fontsize=7.4, labelcolor=INK2,
    handletextpad=0.3, columnspacing=1.0)
letter(axes_c[0], '(c)', x=-0.55, y=1.03)

axd = fig.add_subplot(gs[2, :]); embed(axd, f'{S}/2026-08-20/figures/locus_TNF.png')
letter(axd, '(d)', x=0.10, y=1.02)

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT8}/figures/Fig_adipokine_three_axes.{ext}', dpi=DPI, facecolor=SURFACE)
print('wrote figures/Fig_adipokine_three_axes.png and .pdf')

export_panels(fig, {'a': [axa], 'b': [axb], 'c': axes_c, 'd': [axd]}, 'adipokine_three_axes', f'{OUT8}/figures/panels')
