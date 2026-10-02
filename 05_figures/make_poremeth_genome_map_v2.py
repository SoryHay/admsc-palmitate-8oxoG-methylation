#!/usr/bin/env python3
"""Figure 11a: genome map of the PoreMeth2 segments.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

P = "/data/8oxo_project"
A = f"{P}/AdMSC_DMR"
T = "/data/8oxo_project/AdMSC_paper/sessions/2026-09-03/results/poremeth_map_inputs_COPY_from_job4d9c7639"
OUT = '/data/8oxo_project/AdMSC_paper/sessions/2026-09-03/figures'
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.labelsize": 10, "xtick.labelsize": 9.5, "ytick.labelsize": 10,
    "legend.fontsize": 9.5, "axes.linewidth": 0.8,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

INK, INK2, MUTED = "#1a1a1a", "#4d4d4d", "#8c8c8c"
COL = {"rDNA": "#D55E00", "satellite": "#0072B2", "intergenic": "#009E73", "genic": "#6E6E6E"}
LBL = {"rDNA": "rDNA array", "satellite": "satellite (centromeric / acrocentric)",
       "intergenic": "intergenic, no repeat", "genic": "gene body / promoter"}

length, acc2chr = {}, {}
for line in open(f"{P}/refs/chrom_map.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and f[0] not in ("chrX", "chrY", "chrM"):
        length[f[0]] = int(f[3]); acc2chr[f[1]] = f[0]
CHROMS = [f"chr{i}" for i in range(1, 23)]

cen = {}
for line in open(f"{T}/cen_main.bed"):
    a, st, en = line.split()[:3]
    c = acc2chr.get(a)
    if c:
        cen[c] = (int(st), int(en))

segs = []
for line in open(f"{T}/ann128.tsv"):
    f = line.rstrip("\n").split("\t")
    acc, start, end, meta = f[0], int(f[1]), int(f[2]), f[3]
    gene, prom, isl, rd, ac, ce = (float(x) for x in f[4:10])
    c = acc2chr.get(acc)
    if not c:
        continue
    loc, db, ds, bins, p = meta.split("|")
    db, p = float(db), float(p)
    if rd >= 0.05:
        cls = "rDNA"
    elif ce >= 0.30 or ac >= 0.05:
        cls = "satellite"
    elif gene > 0 or prom > 0:
        cls = "genic"
    else:
        cls = "intergenic"
    segs.append(dict(chrom=c, start=start, end=end, mid=(start + end) / 2,
                     db=db, ds=float(ds), p=p, cls=cls, bins=int(bins), loc=loc))

GENES = {"chr16:92623843-92647536": "FOXC2 / FOXL1",
         "chr19:59916742-60053322": "PEG3 / ZIM2"}

MB = 1e6
xmax = max(length[c] for c in CHROMS) / MB
DBMAX = max(abs(s["db"]) for s in segs) * 1.15
ROW = 1.0

fig, ax = plt.subplots(figsize=(7.4, 8.2))

for i, c in enumerate(CHROMS):
    y = -i * ROW
    L = length[c] / MB

    ax.add_patch(FancyBboxPatch((0, y - 0.07), L, 0.14,
                                boxstyle="round,pad=0,rounding_size=0.35",
                                facecolor="#e8e8e6", edgecolor="none", zorder=1))
    if c in cen:
        cs, ce_ = cen[c][0] / MB, cen[c][1] / MB
        ax.add_patch(FancyBboxPatch((cs, y - 0.07), max(ce_ - cs, 0.4), 0.14,
                                    boxstyle="round,pad=0,rounding_size=0.35",
                                    facecolor="#c2c2be", edgecolor="none", zorder=2))
    ax.text(-2.5, y, c, ha="right", va="center", fontsize=10, color=INK2)

for s in segs:
    if s["db"] == 0:
        continue
    y = -CHROMS.index(s["chrom"]) * ROW
    h = (s["db"] / DBMAX) * 0.42
    sig = s["p"] < 0.05
    ax.plot([s["mid"] / MB, s["mid"] / MB], [y, y + h],
            color=COL[s["cls"]], lw=4.0 if sig else 2.5,
            alpha=0.95 if sig else 0.32, solid_capstyle="round",
            zorder=5 if sig else 3)
    ax.plot([s["mid"] / MB], [y + h], marker="o", ms=8.5 if sig else 5.5,
            mfc=COL[s["cls"]], mec="white", mew=0.9,
            alpha=0.95 if sig else 0.32, zorder=6 if sig else 3)

ax.axhline(0, lw=0)
for i in range(len(CHROMS)):
    ax.plot([0, length[CHROMS[i]] / MB], [-i * ROW, -i * ROW],
            color="#b8b8b4", lw=0.4, zorder=0)

strongest = max(segs, key=lambda s: (-s["p"], abs(s["db"])) if s["p"] < 1e-20 else (-1e9, 0))
ann = []
for s in segs:
    if s["loc"] in GENES:
        ann.append((s, GENES[s["loc"]], -1))
if strongest["p"] < 1e-20:
    ann.append((strongest, f"rDNA, p = {strongest['p']:.0e}", +1))

for s, txt, side in ann:
    y = -CHROMS.index(s["chrom"]) * ROW
    h = (s["db"] / DBMAX) * 0.42
    ax.annotate(txt, xy=(s["mid"] / MB, y + h),
                xytext=(s["mid"] / MB + 14, y + h + (0.30 * side)),
                fontsize=9.5, color=INK,
                arrowprops=dict(arrowstyle="-", lw=0.8, color=MUTED,
                                shrinkA=0, shrinkB=1.5),
                ha="left", va="center", zorder=10)

ax.set_xlim(-14, xmax + 34)
ax.set_ylim(-(len(CHROMS) - 1) * ROW - 0.62, 0.72)
ax.set_yticks([])
ax.set_xlabel("T2T-CHM13v2.0 position (Mb)", color=INK2)
for sp in ("left", "right", "top"):
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color("#b8b8b4")
ax.tick_params(colors=INK2, length=3)
ax.grid(axis="x", color="#ececea", lw=0.5, zorder=-1)
ax.set_axisbelow(True)

handles = [Line2D([], [], color=COL[k], lw=4.0, label=LBL[k])
           for k in ("rDNA", "satellite", "intergenic", "genic")]
handles += [Line2D([], [], color="none", label=""),
            Line2D([], [], color=INK2, lw=1.6, label="up = retains methylation in PA"),
            Line2D([], [], color=INK2, lw=1.6, label="down = loses methylation in PA"),
            Line2D([], [], color=INK2, lw=1.0, alpha=0.32, label="faint = p ≥ 0.05")]
leg = ax.legend(handles=handles, loc="lower right", frameon=False,
                labelspacing=0.38, handlelength=1.6, borderaxespad=0.2)
for t in leg.get_texts():
    t.set_color(INK2)

ax.text(0, 1.030, "Where palmitate changes methylation in AdMSC — PoreMeth2 segments on T2T",
        transform=ax.transAxes, fontsize=8.8, color=INK, va="bottom")
ax.text(0, 1.008, "bar height = Δβ (PA − NT); colour = what the segment overlaps; "
                  "128 segments, 22 autosomes; grey band = centromeric satellite",
        transform=ax.transAxes, fontsize=7.2, color=MUTED, va="bottom")

fig.savefig(f"{OUT}/Fig_poremeth_genome_map_v2.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}/Fig_poremeth_genome_map_v2.pdf", dpi=600, bbox_inches="tight")
print(f"wrote {OUT}/Fig_poremeth_genome_map.png / .pdf")
print(f"segments drawn: {sum(1 for s in segs if s['db'] != 0)} of {len(segs)}")
from collections import Counter
print(Counter(s["cls"] for s in segs))
