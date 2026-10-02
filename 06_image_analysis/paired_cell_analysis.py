"""Cells tracked through all time points; paired per-cell trajectories (Figure 6).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "<DATA_DISK>/imaging_AdMSC_PA/per_cell"
OUT = "/data/imaging_nellie/figures"
ACCENT, INK, GREY = "#1f6f8b", "#22292d", "#b6bfc4"
BASE_T, END_T = (0, 1, 2), 23

df = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{D}/*_stitched.csv"))], ignore_index=True)
n_t = df.t.nunique()

span = df.groupby(["field", "cell_id"]).t.nunique()
complete = set(span[span == n_t].index)
base = df[df.t.isin(BASE_T)]
interior = set(base[~base.touches_border].groupby(["field", "cell_id"]).size().index)
big = set(base.groupby(["field", "cell_id"]).cell_area_um2.mean().pipe(lambda s: s[s > 800]).index)
keep = complete & interior & big
print(f"{len(complete)} complete tracks | interior {len(complete & interior)} | "
      f"and >800 um2: {len(keep)} cells used")
for f in sorted(df.field.unique()):
    print(f"    {f}: {sum(1 for k in keep if k[0]==f)} cells")

df["key"] = list(zip(df.field, df.cell_id))
d = df[df.key.isin(keep)].copy()

PARAMS = [("area_weighted_size_um2", "area-weighted organelle size"),
          ("frac_area_largest", "fraction of area in largest object"),
          ("median_organelle_area_um2", "median organelle area"),
          ("organelle_density_per_100um2", "organelles per 100 µm²"),
          ("mito_area_fraction", "mitochondrial area fraction"),
          ("cv_organelle_area", "CV of organelle size within cell")]
PARAMS = [(c, l) for c, l in PARAMS if c in d.columns]

fig, axes = plt.subplots(2, len(PARAMS), figsize=(3.05 * len(PARAMS), 7.4),
                         gridspec_kw={"height_ratios": [2, 1]})
print("\nPAIRED, each cell against its own 0–2 h baseline:\n")
print(f"{'parameter':28s} {'median':>8s} {'IQR':>14s} {'down/total':>11s} {'sign p':>9s}")
summary = []
for j, (col, label) in enumerate(PARAMS):
    fold = {}
    for k, g in d.groupby("key"):
        b = g[g.t.isin(BASE_T)][col].mean()
        if not np.isfinite(b) or b == 0:
            continue
        fold[k] = (g.set_index("t")[col] / b * 100)
    M = pd.DataFrame(fold)
    ax = axes[0, j]
    for k in M.columns:
        ax.plot(M.index, M[k].values, color=GREY, lw=0.8, alpha=0.75)
    ax.plot(M.index, M.median(axis=1).values, color=ACCENT, lw=2.4)
    ax.axhline(100, color=INK, lw=0.6, ls=(0, (4, 3)), alpha=0.5)
    ax.axvspan(0, 2, color="#c8ced2", alpha=0.3, lw=0)
    ax.set_title(label, fontsize=10, loc="left", color=INK)
    ax.set_xticks([0, 6, 12, 18, 23]); ax.set_ylabel("% of own baseline", fontsize=8, color="#6b757b")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(colors="#6b757b", labelsize=7.5); ax.grid(axis="y", color="#eef1f2", lw=.7)

    end = M.loc[END_T].dropna()
    lo, med, hi = np.percentile(end, [25, 50, 75])
    n_dn = int((end < 100).sum()); n = len(end)
    p = stats.binomtest(max(n_dn, n - n_dn), n, 0.5).pvalue
    summary.append((label, med, lo, hi, n_dn, n, p))
    print(f"{label:28s} {med:7.0f}% {lo:6.0f}–{hi:<6.0f}% {n_dn:5d}/{n:<5d} {p:9.2g}")

    ax2 = axes[1, j]
    ax2.hist(end, bins=np.linspace(min(40, end.min()), max(200, end.max()), 16),
             color=ACCENT, alpha=0.55, edgecolor="white", lw=0.8)
    ax2.axvline(100, color=INK, lw=0.8, ls=(0, (4, 3)))
    ax2.axvline(med, color=ACCENT, lw=2)
    ax2.set_xlabel("% of own baseline at 23 h", fontsize=8, color="#6b757b")
    for s in ("top", "right", "left"): ax2.spines[s].set_visible(False)
    ax2.set_yticks([]); ax2.tick_params(colors="#6b757b", labelsize=7.5)
    dip = stats.shapiro(end).pvalue if 3 <= len(end) <= 5000 else np.nan
    ax2.text(0.98, 0.9, f"n={n}", transform=ax2.transAxes, ha="right", fontsize=8, color=INK)

fig.suptitle("Paired single-cell trajectories — each cell normalised to its own 0–2 h baseline",
             fontsize=13, fontweight="bold", color=INK, x=0.006, ha="left", y=1.0)
fig.text(0.006, 0.965, f"{len(keep)} whole interior cells followed through all 24 timepoints "
         f"by stitched segmentation. Grey, individual cells; blue, median. "
         f"Lower row: distribution of each cell's value at 23 h.",
         fontsize=8.5, color="#6b757b", ha="left", va="top")
fig.tight_layout(rect=(0, 0, 1, 0.945))
os.makedirs(OUT, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fig4_paired_single_cell.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
print(f"\nwrote {OUT}/fig4_paired_single_cell.png")
