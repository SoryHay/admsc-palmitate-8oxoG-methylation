#!/usr/bin/env python3
"""Figure 13.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os, sys
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt

S = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03'
D = f'{S}/p1_5hmC'
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#8a8983'
NT_C, PA_C = '#2a78d6', '#eb6834'

def sums(bed):
    """valid counted once per position (the m and h rows share N_valid_cov); m and h summed."""
    valid = m = h = 0; seen = None
    with open(bed) as f:
        for line in f:
            c = line.split('\t', 12)
            key = (c[0], c[1], c[5])
            if key != seen:
                valid += int(c[9]); seen = key
            if c[3] == 'm': m += int(c[11])
            elif c[3] == 'h': h += int(c[11])
    return valid, m, h

vals = {}
for arm in ('NT', 'PA'):
    v, m, h = sums(f'{D}/{arm}_cpg_thr07.bed')
    vals[arm] = dict(valid=v, m=m, h=h, pm=100 * m / v, ph=100 * h / v)
os.makedirs(f'{S}/results', exist_ok=True)
with open(f'{S}/results/fig_P1_5hmC_v5_values.tsv', 'w') as o:
    o.write('arm\tvalid_calls\t5mC_calls\t5hmC_calls\t5mC_pct\t5hmC_pct\n')
    for arm in ('NT', 'PA'):
        d = vals[arm]; o.write(f"{arm}\t{d['valid']}\t{d['m']}\t{d['h']}\t{d['pm']:.4f}\t{d['ph']:.4f}\n")

from scipy.stats import chi2
fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.9), gridspec_kw={'width_ratios': [1.0, 1.0, 1.25]})
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.07, right=0.99, top=0.90, bottom=0.14, wspace=0.40)
def style(ax):
    ax.set_facecolor('white')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=9, length=3)
    ax.grid(axis='y', color='#ededea', lw=0.8); ax.set_axisbelow(True)

ax = axes[0]; style(ax)
y = [vals['NT']['pm'] / vals['NT']['ph'], vals['PA']['pm'] / vals['PA']['ph']]
ax.plot([0, 1], y, color='#c9c8c3', lw=1.6, zorder=1)
ax.plot([0], [y[0]], 'o', ms=9, color=NT_C, mec='white', mew=1.2, zorder=3)
ax.plot([1], [y[1]], 'o', ms=9, color=PA_C, mec='white', mew=1.2, zorder=3)
top = max(y) * 1.28
for x, v in zip((0, 1), y):
    ax.text(x, v + top * 0.035, f'{v:.1f} : 1', ha='center', va='bottom', fontsize=9, color=INK)
ax.text(0.5, min(y) - top * 0.09, f'×{y[1] / y[0]:.3f}', ha='center', va='top', fontsize=9.5, color=INK2, fontweight='bold')
ax.set_xlim(-0.55, 1.55); ax.set_ylim(0, top)
ax.set_xticks([0, 1]); ax.set_xticklabels(['NT', 'PA'], fontsize=10, color=INK)
ax.set_ylabel('5mC : 5hmC', fontsize=9.5, color=INK2)
ax.text(-0.22, 1.04, '(a)', transform=ax.transAxes, fontsize=12, fontweight='bold', color=INK, ha='left', va='bottom')

ax = axes[1]; style(ax)
d = [vals['PA']['pm'] - vals['NT']['pm'], vals['PA']['ph'] - vals['NT']['ph']]
ax.axhline(0, color='#52514e', lw=0.9, zorder=2)
ax.bar([0, 1], d, width=0.5, color=['#8a8983', '#8a8983'], edgecolor='none', zorder=3)
lim = max(abs(v) for v in d) * 1.6
for x, v in zip((0, 1), d):
    ax.text(x, v + (lim * 0.05 if v > 0 else -lim * 0.05), f'{v:+.2f} pp', ha='center', va='bottom' if v > 0 else 'top', fontsize=9, color=INK)
ax.set_xlim(-0.7, 1.7); ax.set_ylim(-lim, lim)
ax.set_xticks([0, 1]); ax.set_xticklabels(['5mC', '5hmC'], fontsize=10, color=INK)
ax.set_ylabel('change under palmitate, pp of valid CpG calls', fontsize=9, color=INK2)
ax.text(-0.22, 1.04, '(b)', transform=ax.transAxes, fontsize=12, fontweight='bold', color=INK, ha='left', va='bottom')

ax = axes[2]; style(ax)
PTS = [
    (0.82, 'NT', 115.62, 75, 648_674), (1.18, 'PA', 96.29, 80, 830_815),
    (2.32, 'NT', 139.94, 43, 307_257), (2.68, 'PA', 101.34, 41, 404_566)]
for x, arm, rate, k, n in PTS:
    col = NT_C if arm == 'NT' else PA_C
    lo = chi2.ppf(0.025, 2 * k) / 2 / n * 1e6; hi = chi2.ppf(0.975, 2 * (k + 1)) / 2 / n * 1e6
    ax.plot([x, x], [lo, hi], color=col, lw=2, solid_capstyle='round', zorder=2)
    ax.plot(x, rate, 'o', color=col, ms=8, zorder=3)
    ax.text(x + 0.08, rate, f'{rate:.0f}', va='center', fontsize=9, color=INK)
    ax.text(x, lo - 9, f'{k}/{n:,}', ha='center', va='top', fontsize=6.8, color=MUTED)
ax.set_xticks([1.0, 2.5]); ax.set_xticklabels(['ribosomal DNA,\nrRNA class (605 intervals)', 'ribosomal DNA,\nfive acrocentric arrays'], fontsize=9, color=INK)
ax.set_xlim(0.45, 3.05); ax.set_ylim(0, 200)
ax.set_ylabel('8-oxo-dG per million assessable G', fontsize=9.5, color=INK2)
for x, txt in ((1.0, '×0.83  p = 0.26, n.s.'), (2.5, '×0.72  p = 0.15, n.s.')):
    ax.text(x, 6, txt, ha='center', fontsize=8.5, color=INK2)
ax.plot([], [], 'o', color=NT_C, ms=8, label='NT'); ax.plot([], [], 'o', color=PA_C, ms=8, label='PA')
ax.legend(loc='upper right', frameon=False, fontsize=8.5, labelcolor=INK2, handletextpad=0.4)
ax.text(-0.16, 1.04, '(c)', transform=ax.transAxes, fontsize=12, fontweight='bold', color=INK, ha='left', va='bottom')

os.makedirs(f'{S}/figures/for_docx', exist_ok=True)
for ext in ('png', 'pdf'):
    fig.savefig(f'{S}/figures/Fig_P1_5hmC_v5.{ext}', dpi=DPI, facecolor='white')
if DPI >= 400:
    fig.savefig(f'{S}/figures/for_docx/Fig_P1_5hmC_v5.png', dpi=DPI, facecolor='white')
print({a: (round(vals[a]['pm'], 3), round(vals[a]['ph'], 3), vals[a]['valid']) for a in vals})

import sys; sys.path.insert(0, f'{S}/scripts')
from panel_export import export_panels
export_panels(fig, {'a': [axes[0]], 'b': [axes[1]], 'c': [axes[2]]}, '5hmC_and_rDNA_oxidation', f'{S}/figures/panels')
