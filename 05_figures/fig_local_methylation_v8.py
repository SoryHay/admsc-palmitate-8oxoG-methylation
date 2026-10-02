#!/usr/bin/env python3
"""Figure 11.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv, os, math
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-26b_manuscript'
OUT8 = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03'
import os, sys; os.makedirs(f'{OUT8}/figures/for_docx', exist_ok=True); sys.path.insert(0, f'{OUT8}/scripts')
from panel_export import export_panels
R = f'{HERE}/results'
MAP = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03/figures/Fig_poremeth_genome_map_v2.png'
SURFACE, INK, INK2, MUTED = '#ffffff', '#0b0b0b', '#52514e', '#8a8983'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
BG_BODY = 0.8971
FLOOR = 0.10

def rd(f): return list(csv.DictReader(open(f'{R}/{f}'), delimiter='\t'))
def style(ax, xgrid=True):
    ax.set_facecolor(SURFACE)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    ax.grid(axis='x' if xgrid else 'y', color='#ededea', lw=0.8); ax.set_axisbelow(True)

fig = plt.figure(figsize=(13.4, 8.6)); fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[0.95, 1.55],
                      left=0.03, right=0.905, top=0.955, bottom=0.085, wspace=0.30, hspace=0.40)

axa = fig.add_subplot(gs[:, 0]); axa.set_facecolor(SURFACE); axa.axis('off')
img = mpimg.imread(MAP)
h, w = img.shape[:2]
axa.imshow(img[int(h*0.055):int(h*0.985), int(w*0.01):int(w*0.90)], interpolation='lanczos')
axa.text(0, 1.035, '(a)', transform=axa.transAxes, fontsize=13, color=INK, fontweight='semibold', va='bottom')

axb = fig.add_subplot(gs[0, 1]); style(axb)
PANELS = [('adipokine', 'adipokines (13)'), ('oxidative', 'oxidative / redox (64)'),
          ('mito_metabolic', 'mitochondrial / metabolic (83)'), ('senescence', 'ageing / senescence (75)'),
          ('innate_exons', 'innate immunity (161)')]
for y, (key, label) in enumerate(PANELS):
    rows = rd(f'layer1_cpg_{key}_bygene_COPY_from_AdMSC_DMR.tsv')
    vals = [float(r['ratio_PA_NT']) / BG_BODY for r in rows
            if float(r['NT_rate']) >= FLOOR and float(r['ratio_PA_NT']) > 0]
    greyed = sum(1 for r in rows if float(r['NT_rate']) < FLOOR)
    xs = [math.log2(v) for v in vals]
    import random; random.seed(11)
    axb.scatter(xs, [y + random.uniform(-0.22, 0.22) for _ in xs], s=13, color=BLUE, alpha=0.55,
                edgecolors='none', zorder=3)
    med = sorted(xs)[len(xs) // 2]
    axb.plot([med, med], [y - 0.32, y + 0.32], color=INK, lw=2, zorder=4)
axb.axvline(0, color='#6f6e69', lw=1.3, zorder=2)
axb.set_yticks(range(len(PANELS))); axb.set_yticklabels([p[1] for p in PANELS], fontsize=8.5)
axb.set_ylim(-0.6, len(PANELS) - 0.4); axb.set_xlim(-2.1, 2.1)
axb.set_xticks([-2, -1, 0, 1, 2]); axb.set_xticklabels(['0.25', '0.5', '1', '2', '4'])
axb.set_xlabel('per-gene ratio relative to the gene-body background 0.8971 (log scale)',
               fontsize=9, color=INK)
axb.text(-0.26, 1.02, '(b)', transform=axb.transAxes, fontsize=13, color=INK, fontweight='semibold', va='bottom')

axc = fig.add_subplot(gs[1, 1]); style(axc)
ax_rows = rd('d11_inflammasome_axis_COPY_from_2026-08-14.tsv')
GROUP = {
    'substrate': ('IL-1β production', BLUE), 'protease': ('IL-1β production', BLUE),
    'extracellular processing': ('IL-1β production', BLUE), 'co-substrate': ('IL-1β production', BLUE),
    'sensor': ('sensing and release', AQUA), 'pore': ('sensing and release', AQUA),
    'adaptor': ('sensing and release', AQUA), 'licensing': ('sensing and release', AQUA),
    'priming': ('priming', ORANGE),
    'non-canonical': ('other', YELLOW), 'alternative': ('other', YELLOW), 'receptor': ('other', YELLOW),
}
def group(role):
    for k, v in GROUP.items():
        if k in role: return v
    return ('other', YELLOW)
order_groups = ['IL-1β production', 'priming', 'sensing and release', 'other']
rows = []
for r in ax_rows:
    if r.get('scale', 'span') != 'span' or not r['rel_to_matched_background']:
        continue
    g, c = group(r['role'])
    rows.append((order_groups.index(g), float(r['rel_to_matched_background']), r, g, c))
rows.sort(key=lambda t: (t[0], -t[1]))
yticks, ylabels, seen = [], [], set()
for y, (gi, rel, r, g, c) in enumerate(rows):
    nt = float(r['NT_pct']); calls = min(int(r['NT_calls']), int(r['PA_calls']))
    uninterp = nt < 10.0
    x = max(rel, 0.05)
    axc.plot(x, y, 'o', ms=7.5 if calls >= 100 else 5.5, zorder=4,
             color=SURFACE if uninterp else c, markeredgecolor=c, markeredgewidth=1.6)
    lab = r['gene'] + ('  ·' if calls < 100 else '')
    yticks.append(y); ylabels.append(lab)
    if g not in seen:
        seen.add(g)
        axc.text(1.012, 1 - (y + 0.5) / len(rows), g, transform=axc.transAxes, va='center', ha='left',
                 fontsize=8.5, color=c, fontweight='medium')
axc.axvline(1.0, color='#6f6e69', lw=1.3, zorder=2)
axc.set_yticks(yticks); axc.set_yticklabels(ylabels, fontsize=7.6)
axc.invert_yaxis(); axc.set_xlim(0.2, 2.3); axc.set_xscale('log')
axc.set_xticks([0.25, 0.5, 1, 2]); axc.set_xticklabels(['0.25', '0.5', '1', '2'])
axc.set_xlabel('gene-body ratio relative to its matched background (log scale)\n'
               '· after a name: fewer than 100 calls   ·   open marker: NT methylation below 10 %',
               fontsize=8.6, color=INK)
axc.text(-0.26, 1.02, '(c)', transform=axc.transAxes, fontsize=13, color=INK, fontweight='semibold', va='bottom')

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT8}/figures/Fig_local_methylation.{ext}', dpi=DPI, facecolor=SURFACE)
print('wrote figures/Fig_local_methylation.png and .pdf')

export_panels(fig, {'a': [axa], 'b': [axb], 'c': [axc]}, 'local_methylation', f'{OUT8}/figures/panels')
