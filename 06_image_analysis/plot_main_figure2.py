"""Figure 5.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations
import glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

DATA = "<DATA_DISK>/imaging_AdMSC_PA"
OUT = "/data/imaging_nellie/figures_manuscript"
FIELDS = ("F0", "F1", "F2")
BASE_T = [0, 1, 2]
ACCENT, FIELD_GREY, NOISE_GREY, INK = "#1f6f8b", "#9aa5ab", "#c8ced2", "#22292d"

PANELS = [("organelle_solidity_mean", "organelle solidity", "compactness"),
          ("organelle_extent_mean", "organelle extent", "compactness"),
          ("branch_aspect_ratio_mean", "branch aspect ratio", "elongation"),
          ("branch_tortuosity_mean", "branch tortuosity", "elongation"),
          ("branch_solidity_mean", "branch solidity", "compactness"),
          ("organelle_area_mean", "organelle area", "size"),
          ("branch_length_mean", "branch length", "length"),
          ("branch_thickness_mean", "branch thickness", "calibre — unchanged"),
          ("intensity_mean", "intensity", "signal — not interpreted")]

d = {f: pd.read_csv(glob.glob(f"{DATA}/prod_deconvolved_{f}/*features_image.csv")[0]).set_index("t").sort_index()
     for f in FIELDS}
hours = np.arange(24)

fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.0))
for ax, (col, label, family) in zip(axes.ravel(), PANELS):
    norm, basepts = [], []
    for f in FIELDS:
        v = d[f][col].astype(float); b = v.loc[BASE_T].mean()
        r = v / b * 100; norm.append(r); basepts += list(r.loc[BASE_T])
    M = pd.concat(norm, axis=1); mean, sd = M.mean(axis=1), M.std(axis=1)
    floor = float(np.std(basepts, ddof=1))

    ax.axhspan(100 - floor, 100 + floor, color=NOISE_GREY, alpha=0.6, lw=0, zorder=0)
    ax.axvspan(0, 2, color=NOISE_GREY, alpha=0.3, lw=0, zorder=0)
    ax.axhline(100, color=INK, lw=0.6, ls=(0, (4, 3)), alpha=0.5, zorder=1)
    for r in norm:
        ax.plot(hours, r.values, color=FIELD_GREY, lw=1.0, alpha=0.85, zorder=2)
    ax.fill_between(hours, (mean - sd).values, (mean + sd).values, color=ACCENT, alpha=0.20, lw=0, zorder=3)
    ax.plot(hours, mean.values, color=ACCENT, lw=2.3, zorder=4, solid_capstyle="round")
    end = mean.iloc[-1]
    ax.annotate(f"{end:.0f}%", xy=(23, end), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=10, color=INK, fontweight="bold")
    ax.set_title(label, fontsize=11.5, loc="left", color=INK, pad=6)

    null = "unchanged" in family or "not interpreted" in family
    ax.text(0.02, 0.055, f"{family}   ·   floor ±{floor:.0f}%   ·   "
            f"{abs(end-100)/floor:.1f}× floor",
            transform=ax.transAxes, fontsize=8.2,
            color="#b5651d" if null else "#6b757b",
            fontweight="bold" if null else "normal")
    ax.set_xlim(-0.6, 26); ax.set_xticks([0, 6, 12, 18, 23])
    ax.grid(axis="y", color="#e6eaec", lw=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color("#c3cacd")
    ax.tick_params(colors="#6b757b", labelsize=8.5)
for ax in axes[-1]: ax.set_xlabel("hours after palmitate", fontsize=9, color="#6b757b")
for ax in axes[:, 0]: ax.set_ylabel("% of baseline", fontsize=9, color="#6b757b")

handles = [plt.Line2D([], [], color=ACCENT, lw=2.3, label="mean of 3 fields"),
           plt.Rectangle((0, 0), 1, 1, color=ACCENT, alpha=0.20, label="± SD between fields"),
           plt.Line2D([], [], color=FIELD_GREY, lw=1.0, label="individual fields"),
           plt.Rectangle((0, 0), 1, 1, color=NOISE_GREY, alpha=0.6, label="technical floor"),
           plt.Rectangle((0, 0), 1, 1, color=NOISE_GREY, alpha=0.3, label="0–2 h pre-response")]
fig.legend(handles=handles, loc="upper left", ncol=5, frameon=False, fontsize=9,
           bbox_to_anchor=(0.007, 0.995))
fig.suptitle("Mitochondrial network morphology over 24 h of palmitate",
             fontsize=13.5, fontweight="bold", color=INK, x=0.007, ha="left", y=1.045)
fig.tight_layout(rect=(0, 0, 1, 0.945))
os.makedirs(OUT, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/Figure2_network_morphology.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
print("wrote Figure2_network_morphology.png/.pdf (9 panels)")
