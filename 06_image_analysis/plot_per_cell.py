"""Figure S3.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

D = "<DATA_DISK>/imaging_AdMSC_PA/per_cell"
OUT = "/data/imaging_nellie/figures"
ACCENT, INK, GREY = "#1f6f8b", "#22292d", "#9aa5ab"

df = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{D}/*_per_cell.csv"))], ignore_index=True)
print(f"{len(df)} cell-timepoint rows, {df.field.nunique()} fields, "
      f"{df.groupby(['field','t']).size().mean():.0f} cells per field-timepoint on average")

PARAMS = [
    ("mito_area_fraction", "mitochondrial area fraction of cell", True),
    ("mean_organelle_area_um2", "mean organelle area (µm²)", False),
    ("organelle_density_per_100um2", "organelles per 100 µm² of cell", True),
    ("branch_length_mean", "mean branch length (µm)", False),
    ("organelle_solidity_raw", "organelle solidity", False),
    ("branch_solidity_mean", "branch solidity", False),
    ("branch_aspect_ratio_mean", "branch aspect ratio", False),
]
PARAMS = [(c, l, b) for c, l, b in PARAMS if c in df.columns]

fig, axes = plt.subplots(2, 3, figsize=(16, 8.6))
ts = sorted(df.t.unique())
for ax, (col, label, drop_border) in zip(axes.ravel(), PARAMS):
    sub = df[~df.touches_border] if drop_border else df
    data = [sub.loc[sub.t == t, col].dropna().values for t in ts]
    parts = ax.violinplot(data, positions=ts, widths=0.85, showextrema=False, showmedians=False)
    for b in parts["bodies"]:
        b.set_facecolor(ACCENT); b.set_alpha(0.30); b.set_edgecolor("none")
    med = [np.median(d) if len(d) else np.nan for d in data]
    q1 = [np.percentile(d, 25) if len(d) else np.nan for d in data]
    q3 = [np.percentile(d, 75) if len(d) else np.nan for d in data]
    ax.fill_between(ts, q1, q3, color=ACCENT, alpha=0.22, lw=0)
    ax.plot(ts, med, color=ACCENT, lw=2.2, zorder=5)
    ax.axvspan(0, 2, color="#c8ced2", alpha=0.30, lw=0, zorder=0)
    ax.set_title(label + ("" if drop_border else "  (all cells)"), fontsize=10, loc="left", color=INK)
    ax.set_xticks([0, 6, 12, 18, 23]); ax.set_xlabel("hours after palmitate", fontsize=8, color="#6b757b")
    ax.grid(axis="y", color="#e6eaec", lw=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(colors="#6b757b", labelsize=8)
    b0 = np.nanmedian(np.concatenate([data[t] for t in (0, 1, 2) if len(data[t])]))
    ax.annotate(f"{med[-1]/b0*100:.0f}% of 0–2 h median", xy=(0.98, 0.94),
                xycoords="axes fraction", ha="right", fontsize=8, color=INK, fontweight="bold")

fig.suptitle("Per-cell distributions across 24 h of palmitate", fontsize=13.5, x=0.008,
             ha="left", y=0.985, fontweight="bold", color=INK)
fig.text(0.008, 0.955, "Each violin is the population of cells at that hour (3 fields pooled); "
         "line = median, band = interquartile range. Cells are not tracked between timepoints.",
         fontsize=8.6, color="#6b757b", ha="left", va="top")
fig.tight_layout(rect=(0, 0, 1, 0.93))
os.makedirs(OUT, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fig3_per_cell_distributions.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
print(f"wrote {OUT}/fig3_per_cell_distributions.png")

print("\nuniform shift or subpopulation? (robust CV = IQR/median, per timepoint, pooled fields)")
for col, label, _ in PARAMS:
    s = df[~df.touches_border] if col in ("mito_area_fraction", "organelle_density_per_100um2") else df
    g = s.groupby("t")[col]
    cv = ((g.quantile(.75) - g.quantile(.25)) / g.median())
    print(f"  {label:38s} 0-2 h {cv.loc[[0,1,2]].mean():.3f}  ->  23 h {cv.loc[23]:.3f}"
          f"   ({cv.loc[23]/cv.loc[[0,1,2]].mean():.2f}x)")
