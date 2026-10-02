#!/usr/bin/env python3
"""Figure 8.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv
import os
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.stats import chi2

HERE = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-26b_manuscript'
OUT8 = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03'
import os, sys; os.makedirs(f'{OUT8}/figures/for_docx', exist_ok=True); sys.path.insert(0, f'{OUT8}/scripts')
from panel_export import export_panels
R = f'{HERE}/results'
SURFACE, INK, INK2, MUTED = '#ffffff', '#0b0b0b', '#52514e', '#8a8983'
NT_C, PA_C = '#2a78d6', '#eb6834'
BIN = 500

def rd(f):
    return list(csv.DictReader(open(f'{R}/{f}'), delimiter='\t'))

def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    ax.grid(axis='y', color='#ededea', lw=0.8)
    ax.set_axisbelow(True)

fig = plt.figure(figsize=(11.6, 10.6))
fig.patch.set_facecolor(SURFACE)
outer = fig.add_gridspec(2, 1, height_ratios=[1.00, 2.25], hspace=0.31,
                         left=0.095, right=0.975, top=0.965, bottom=0.062)

ax = fig.add_subplot(outer[0]); style(ax)
PTS = [
    (0.82, 'NT', 95.45, 34, 356_221),  (1.18, 'PA', 121.15, 28, 231_114),
    (2.32, 'NT', 62.10, 48, 772_996),  (2.68, 'PA', 49.37, 45, 911_550)]
for x, arm, rate, k, n in PTS:
    col = NT_C if arm == 'NT' else PA_C
    lo = chi2.ppf(0.025, 2 * k) / 2 / n * 1e6
    hi = chi2.ppf(0.975, 2 * (k + 1)) / 2 / n * 1e6
    ax.plot([x, x], [lo, hi], color=col, lw=2, solid_capstyle='round', zorder=2)
    ax.plot(x, rate, 'o', color=col, ms=9, zorder=3)
    ax.text(x + 0.08, rate, f'{rate:.0f}', va='center', fontsize=10.5, color=INK)
    ax.text(x, lo - 9, f'{k}/{n:,}', ha='center', va='top', fontsize=7.6, color=MUTED)
ax.set_xticks([1.0, 2.5])
ax.set_xticklabels(['mitochondrial genome', 'nuclear 82-gene panel'], fontsize=10.5, color=INK)
ax.set_xlim(0.45, 3.05); ax.set_ylim(0, 200)
ax.set_ylabel('8-oxo-dG per million\nassessable guanines', fontsize=9.5, color=INK)
for x, txt in ((1.0, 'PA / NT = 1.27×   p = 0.36, n.s.'), (2.5, 'PA / NT = 0.80×   p = 0.30, n.s.')):
    ax.text(x, 8, txt, ha='center', fontsize=9.5, color=INK2)
ax.plot([], [], 'o', color=NT_C, ms=9, label='untreated (NT)')
ax.plot([], [], 'o', color=PA_C, ms=9, label='palmitate (PA)')
ax.legend(loc='upper right', frameon=False, fontsize=9.5, labelcolor=INK2, handletextpad=0.4)

inner = outer[1].subgridspec(3, 1, height_ratios=[0.30, 1.0, 1.15], hspace=0.10)
axg = fig.add_subplot(inner[0])
ax1 = fig.add_subplot(inner[1], sharex=axg)
ax2 = fig.add_subplot(inner[2], sharex=axg)
for a in (ax1, ax2):
    style(a)
axg.set_facecolor(SURFACE); axg.axis('off')

GENES = [('RNR1', 72, 1025), ('RNR2', 1095, 2652), ('ND1', 2730, 3685), ('ND2', 3893, 4934),
         ('CO1', 5327, 6868), ('CO2', 7009, 7692), ('ATP6', 7950, 8630), ('CO3', 8630, 9413),
         ('ND4', 10183, 11560), ('ND5', 11760, 13571), ('ND6', 13572, 14096), ('CYB', 14170, 15310)]
axg.add_patch(Rectangle((0, 0.34), 16569, 0.32, facecolor='#eceae4', edgecolor='none'))
for name, a0, a1 in GENES:
    axg.add_patch(Rectangle((a0, 0.22), a1 - a0, 0.56, facecolor='#a8a69c', edgecolor='none'))
    if a1 - a0 > 600:
        axg.text((a0 + a1) / 2, 0.50, name, ha='center', va='center', fontsize=7.5, color='white')
axg.set_xlim(0, 16569); axg.set_ylim(0, 1)
axg.text(0, 1.55, 'mitochondrial genome, T2T CP068254.1 (16,569 bp)', fontsize=8.8, color=MUTED)

rows = sorted(rd('d6_mito_bins_COPY_from_2026-08-14.tsv'), key=lambda r: int(r['bin_start']))
agg = {}
for r in rows:
    b = (int(r['bin_start']) - 1) // BIN
    d = agg.setdefault(b, dict(cg_NT=0, c_NT=0, cg_PA=0, c_PA=0))
    for arm in ('NT', 'PA'):
        d[f'cg_{arm}'] += int(r[f'callable_G_{arm}'])
        d[f'c_{arm}'] += int(r[f'calls_0.70_{arm}'])
xs = [(b + 0.5) * BIN for b in sorted(agg)]
ynt = [agg[b]['cg_NT'] for b in sorted(agg)]; ypa = [agg[b]['cg_PA'] for b in sorted(agg)]
ax1.fill_between(xs, ypa, ynt, color='#e3e2dc', alpha=0.9, lw=0, zorder=1)
for arm, col, ys in (('NT', NT_C, ynt), ('PA', PA_C, ypa)):
    ax1.plot(xs, ys, color=col, lw=1.8, zorder=3)
    ax2.plot(xs, [1e5 * agg[b][f'c_{arm}'] / max(1, agg[b][f'cg_{arm}']) for b in sorted(agg)],
             color=col, lw=1.8, marker='o', ms=2.6, mew=0)
ax1.set_ylabel('assessable G\nper 500 bp', fontsize=9.5, color=INK)
ax2.set_ylabel('8-oxo-dG calls per 100,000\nassessable G (relaxed cut 0.70)', fontsize=9.5, color=INK)
ax1.tick_params(labelbottom=False); axg.tick_params(labelbottom=False)

calls = rd('d6_mito_calls_locked_COPY_from_2026-08-14.tsv')
n_nt = sum(1 for c in calls if c['arm'] == 'NT')
n_pa = len(calls) - n_nt
for arm, col, y in (('NT', NT_C, -95), ('PA', PA_C, -155)):
    for c in calls:
        if c['arm'] == arm:
            ax2.plot([int(c['chrM_pos'])] * 2, [y - 26, y + 26], color=col, lw=1.6)
ax2.text(16300, -95, f'NT ({n_nt})', color=NT_C, fontsize=8, va='center', ha='right')
ax2.text(16300, -155, f'PA ({n_pa})', color=PA_C, fontsize=8, va='center', ha='right')
ax2.set_ylim(-200, 1080)
ax2.set_xlabel('position on the mitochondrial genome (bp)', fontsize=9.5, color=INK)
ax2.set_xlim(0, 16569)
ax2.text(200, -70, 'ticks below: every call at the locked per-k-mer operating point',
         fontsize=8, color=MUTED, va='bottom')

fig.text(0.012, 0.655, '(b)', fontsize=13, color=INK, fontweight='semibold', va='top')
fig.text(0.012, 0.958, '(a)', fontsize=13, color=INK, fontweight='semibold', va='top')

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT8}/figures/Fig_mito_oxidation_combined.{ext}', dpi=DPI, facecolor=SURFACE)
print('wrote figures/Fig_mito_oxidation_combined.png and .pdf')

export_panels(fig, {'a': [ax], 'b': [axg, ax1, ax2]}, 'mito_oxidation', f'{OUT8}/figures/panels')
