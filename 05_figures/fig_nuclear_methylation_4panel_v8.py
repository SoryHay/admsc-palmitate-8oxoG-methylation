#!/usr/bin/env python3
"""Figure 9.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv
import os
import matplotlib
matplotlib.use('Agg')
DPI = int(os.environ.get('FIG_DPI', 200))
import matplotlib.pyplot as plt

HERE = '/data/8oxo_project/AdMSC_paper/sessions/2026-08-26b_manuscript'
OUT8 = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03'
import os, sys; os.makedirs(f'{OUT8}/figures/for_docx', exist_ok=True); sys.path.insert(0, f'{OUT8}/scripts')
from panel_export import export_panels
R = f'{HERE}/results'
SURFACE, INK, INK2, MUTED = '#ffffff', '#0b0b0b', '#52514e', '#8a8983'
NT_C, PA_C, HL = '#2a78d6', '#eb6834', '#1baf7a'
GENOME = 0.89845
FLOOR = 10.0

def rd(f):
    return list(csv.DictReader(open(f'{R}/{f}'), delimiter='\t'))

def style(ax, xgrid=False):
    ax.set_facecolor(SURFACE)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#d8d7d2')
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    ax.grid(axis='x' if xgrid else 'y', color='#ededea', lw=0.8)
    ax.set_axisbelow(True)

def letter(ax, s, x=-0.115, y=1.10):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=13, color=INK,
            fontweight='semibold', va='top', ha='left')

fig = plt.figure(figsize=(12.6, 9.4))
fig.patch.set_facecolor(SURFACE)
outer = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30,
                         left=0.085, right=0.975, top=0.955, bottom=0.075)

ax = fig.add_subplot(outer[0, 0]); style(ax, xgrid=True)
rows = sorted(rd('f1_global_shift_COPY_from_2026-08-14.tsv'), key=lambda r: float(r['NT_pct']))
ys = range(len(rows))
for y, r in zip(ys, rows):
    nt, pa = float(r['NT_pct']), float(r['PA_pct'])
    ax.plot([pa, nt], [y, y], color='#cdccc6', lw=2.2, zorder=1, solid_capstyle='round')
    ax.plot(nt, y, 'o', color=NT_C, ms=7, zorder=3)
    ax.plot(pa, y, 'o', color=PA_C, ms=7, zorder=3)
    ax.text(45.5, y, f"{float(r['ratio_PA_NT']):.4f}", va='center', ha='right',
            fontsize=8.5, color=INK)
ax.set_yticks(list(ys))
ax.set_yticklabels([r['region_class'].replace(' (TSS+-1kb)', '') for r in rows], fontsize=8.5)
for t, r in zip(ax.get_yticklabels(), rows):
    if r['region_class'] == 'whole genome':
        t.set_fontweight('bold'); t.set_color(INK)
ax.set_xlim(0, 46); ax.set_ylim(-0.7, len(rows) - 0.3)
ax.set_xlabel('CpG methylation (%)', fontsize=9.5, color=INK)
ax.text(45.5, len(rows) - 0.35, 'PA / NT', ha='right', fontsize=8.5, color=INK2, fontweight='medium')
ax.plot([], [], 'o', color=NT_C, ms=7, label='untreated (NT)')
ax.plot([], [], 'o', color=PA_C, ms=7, label='palmitate (PA)')
ax.legend(loc='upper left', frameon=False, fontsize=8.5, labelcolor=INK2,
          handletextpad=0.3, borderpad=0.1)
letter(ax, '(a)'); ax_a = ax

ax = fig.add_subplot(outer[0, 1]); style(ax)
b = rd('f2_per_read_beta_COPY_from_2026-08-14.tsv')
for arm, col, mk in (('NT', NT_C, 'o'), ('PA', PA_C, 's')):
    s = sorted((r for r in b if r['arm'] == arm), key=lambda r: float(r['beta_lo']))
    x = [min(float(r['beta_lo']), 1.0) + 0.025 for r in s]
    y = [100 * float(r['fraction_of_reads']) for r in s]
    ax.plot(x, y, color=col, lw=2, marker=mk, ms=4, mew=0, label=f'{arm}')
ax.set_xlim(-0.03, 1.09); ax.set_ylim(0, 32)
ax.set_xlabel('methylation of the read (β)', fontsize=9.5, color=INK)
ax.set_ylabel('share of reads (%)', fontsize=9.5, color=INK)
ax.legend(loc='upper right', frameon=False, fontsize=8.5, labelcolor=INK2, handletextpad=0.4)
ax.annotate('unmethylated reads\n26.1 → 29.2 %', xy=(0.025, 29), xytext=(0.15, 25.5),
            fontsize=8.5, color=INK2,
            arrowprops=dict(arrowstyle='->', color=INK2, lw=1))
ax.annotate('reads above β 0.6\n14.8 → 12.2 %', xy=(0.72, 2.5), xytext=(0.52, 9),
            fontsize=8.5, color=INK2,
            arrowprops=dict(arrowstyle='->', color=INK2, lw=1))
letter(ax, '(b)'); ax_b = ax

ax = fig.add_subplot(outer[1, 0]); style(ax, xgrid=True)
rr = sorted(rd('f6_repeat_classes_COPY_from_2026-08-14.tsv'), key=lambda r: float(r['ratio_PA_NT']))
ax.axvline(1.0, color='#6f6e69', lw=1.3, zorder=2)
ax.axvline(GENOME, color='#9c9b95', lw=1.2, ls=(0, (5, 4)), zorder=2)
for y, r in enumerate(rr):
    v, ntp = float(r['ratio_PA_NT']), float(r['NT_pct'])
    is_rdna = r['repeat_class'] == 'rRNA'
    col = HL if is_rdna else NT_C
    ax.plot(v, y, 'o', ms=8 if is_rdna else 6.5, zorder=4,
            color=col if ntp >= FLOOR else SURFACE,
            markeredgecolor=col, markeredgewidth=1.8)
ax.set_yticks(range(len(rr)))
ax.set_yticklabels([('rDNA' if r['repeat_class'] == 'rRNA' else r['repeat_class'].replace('_', ' '))
                    for r in rr], fontsize=8.5)
for t, r in zip(ax.get_yticklabels(), rr):
    if r['repeat_class'] == 'rRNA':
        t.set_color('#12805a'); t.set_fontweight('bold')
ax.set_xlim(0.855, 1.045); ax.set_ylim(-0.7, len(rr) - 0.3)
ax.set_xlabel('PA / NT methylation ratio        (open marker: NT methylation below 10 %)',
              fontsize=9.5, color=INK)
ax.text(GENOME - 0.0025, len(rr) - 0.6, 'whole genome 0.898', rotation=90, ha='right', va='top',
        fontsize=8, color=INK2)
ax.text(1.0025, len(rr) - 0.6, 'no change', rotation=90, ha='left', va='top',
        fontsize=8, color=INK2)
letter(ax, '(c)'); ax_c = ax

inner = outer[1, 1].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.10)
ax1 = fig.add_subplot(inner[0]); ax2 = fig.add_subplot(inner[1], sharex=ax1)
for a in (ax1, ax2):
    style(a)
m = sorted(rd('f7_metagene_COPY_from_2026-08-14.tsv'), key=lambda r: int(r['bin']))
x = [int(r['bin']) for r in m]
ax1.plot(x, [float(r['NT_pct']) for r in m], color=NT_C, lw=2, label='NT')
ax1.plot(x, [float(r['PA_pct']) for r in m], color=PA_C, lw=2, label='PA')
ax2.plot(x, [float(r['ratio']) for r in m], color=INK, lw=1.4)
ax2.axhline(GENOME, color='#9c9b95', lw=1.2, ls=(0, (5, 4)))
for a in (ax1, ax2):
    for xb in (19.5, 119.5):
        a.axvline(xb, color='#cdccc6', lw=1)
    a.set_xlim(0, 139)
ax1.set_ylim(0, 34); ax2.set_ylim(0.855, 0.945)
ax1.tick_params(labelbottom=False)
ax1.set_ylabel('CpG methylation (%)', fontsize=9.5, color=INK)
ax2.set_ylabel('PA / NT', fontsize=9.5, color=INK)
ax2.set_xticks([0, 19.5, 69, 119.5, 139])
ax2.set_xticklabels(['−2 kb', 'TSS', 'gene body', 'TES', '+2 kb'], fontsize=8.5)
ax1.legend(loc='lower right', frameon=False, fontsize=8.5, labelcolor=INK2, handletextpad=0.4)
ax2.text(2, 0.861, 'dashed: whole genome 0.898', ha='left', fontsize=8, color=INK2)
letter(ax1, '(d)', y=1.16)

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT8}/figures/Fig_nuclear_methylation_4panel.{ext}', dpi=DPI, facecolor=SURFACE)
print('wrote figures/Fig_nuclear_methylation_4panel.png and .pdf')

export_panels(fig, {'a': [ax_a], 'b': [ax_b], 'c': [ax_c], 'd': [ax1, ax2]}, 'nuclear_methylation', f'{OUT8}/figures/panels')
