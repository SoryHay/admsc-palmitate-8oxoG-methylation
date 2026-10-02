"""Figures S1-S2.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = "<DATA_DISK>/imaging_AdMSC_PA"
OUT = "/data/imaging_nellie/figures"
FIELDS = ("F0", "F1", "F2")
BASELINE_T = (0, 1, 2)

ACCENT = "#1f6f8b"
FIELD_GREY = "#9aa5ab"
NOISE_GREY = "#c8ced2"
INK = "#22292d"

RATIO_SUMS = {
    "branch_tortuosity_sum", "branch_solidity_sum", "branch_extent_sum",
    "branch_aspect_ratio_sum", "organelle_solidity_sum", "organelle_extent_sum",
    "branch_thickness_sum", "branch_axis_length_maj_sum", "branch_axis_length_min_sum",
    "organelle_axis_length_maj_sum", "organelle_axis_length_min_sum",
}

PRETTY = {
    "structure": "structure", "intensity": "intensity",
    "branch_length": "branch length", "branch_thickness": "branch thickness",
    "branch_aspect_ratio": "branch aspect ratio", "branch_tortuosity": "branch tortuosity",
    "branch_area": "branch area", "branch_axis_length_maj": "branch major axis",
    "branch_axis_length_min": "branch minor axis", "branch_extent": "branch extent",
    "branch_solidity": "branch solidity", "organelle_area": "organelle area",
    "organelle_axis_length_maj": "organelle major axis",
    "organelle_axis_length_min": "organelle minor axis",
    "organelle_extent": "organelle extent", "organelle_solidity": "organelle solidity",
}

def load():
    d = {}
    for f in FIELDS:
        p = glob.glob(f"{DATA}/prod_deconvolved_{f}/*features_image.csv")[0]
        d[f] = pd.read_csv(p).set_index("t").sort_index()
    return d

def panel(ax, data, col, hours):
    norm, base_pts = [], []
    for f in FIELDS:
        v = data[f][col].astype(float)
        base = v.loc[list(BASELINE_T)].mean()
        if not np.isfinite(base) or base == 0:
            ax.text(0.5, 0.5, "no baseline", ha="center", va="center",
                    transform=ax.transAxes, color=FIELD_GREY, fontsize=7)
            return None
        r = v / base * 100.0
        norm.append(r)
        base_pts.extend(r.loc[list(BASELINE_T)].tolist())

    M = pd.concat(norm, axis=1)
    mean, sd = M.mean(axis=1), M.std(axis=1)

    noise = float(np.std(base_pts, ddof=1))
    ax.axhspan(100 - noise, 100 + noise, color=NOISE_GREY, alpha=0.55, lw=0, zorder=0)
    ax.axvspan(0, 2, color=NOISE_GREY, alpha=0.30, lw=0, zorder=0)
    ax.axhline(100, color=INK, lw=0.6, ls=(0, (4, 3)), alpha=0.5, zorder=1)

    for r in norm:
        ax.plot(hours, r.values, color=FIELD_GREY, lw=0.9, alpha=0.85, zorder=2)
    ax.fill_between(hours, (mean - sd).values, (mean + sd).values,
                    color=ACCENT, alpha=0.20, lw=0, zorder=3)
    ax.plot(hours, mean.values, color=ACCENT, lw=2.0, zorder=4,
            solid_capstyle="round")

    end = mean.iloc[-1]
    ax.annotate(f"{end:.0f}%", xy=(hours[-1], end), xytext=(3, 0),
                textcoords="offset points", va="center", ha="left",
                fontsize=7.5, color=INK, fontweight="bold")
    return noise

def make_figure(data, cols, title, subtitle, outstem, mark_ratio=False):
    hours = np.arange(24)
    n = len(cols)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.5, 2.9 * nrow), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, col in zip(axes, cols):
        noise = panel(ax, data, col, hours)
        stem = col.rsplit("_", 1)[0]
        name = PRETTY.get(stem, stem.replace("_", " "))
        flag = "  †" if (mark_ratio and col in RATIO_SUMS) else ""
        ax.set_title(f"{name}{flag}", fontsize=9.5, color=INK, pad=6, loc="left")
        if noise is not None:

            ax.text(0.02, 0.04, f"floor ±{noise:.0f}%", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=6.5, color="#6b757b")
        ax.set_xlim(-0.6, 25.5)
        ax.set_xticks([0, 6, 12, 18, 23])
        ax.grid(axis="y", color="#e6eaec", lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#c3cacd")
        ax.tick_params(colors="#6b757b", labelsize=7.5, length=3)

    for ax in axes[n:]:
        ax.set_visible(False)
    for ax in axes[max(0, n - ncol):n]:
        ax.set_xlabel("hours after palmitate", fontsize=8, color="#6b757b")
    for i in range(0, n, ncol):
        axes[i].set_ylabel("% of baseline", fontsize=8, color="#6b757b")

    head = 1.0 - 0.30 / (2.9 * nrow)
    fig.suptitle(title, fontsize=13.5, color=INK, x=0.008, ha="left", y=head + 0.012,
                 fontweight="bold")
    fig.text(0.008, head - 0.004, subtitle, fontsize=8.6, color="#6b757b", ha="left",
             va="top", wrap=True)

    handles = [
        plt.Line2D([], [], color=ACCENT, lw=2.0, label="mean of 3 fields"),
        plt.Rectangle((0, 0), 1, 1, color=ACCENT, alpha=0.20, label="± SD between fields"),
        plt.Line2D([], [], color=FIELD_GREY, lw=0.9, label="individual fields (F0, F1, F2)"),
        plt.Rectangle((0, 0), 1, 1, color=NOISE_GREY, alpha=0.55,
                      label="baseline technical floor (± SD of 0–2 h)"),
        plt.Rectangle((0, 0), 1, 1, color=NOISE_GREY, alpha=0.30,
                      label="0–2 h pre-response window"),
    ]
    fig.legend(handles=handles, loc="upper left", ncol=5, frameon=False,
               fontsize=8, bbox_to_anchor=(0.008, head - 0.030))

    fig.tight_layout(rect=(0, 0, 1, head - 0.048))
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{outstem}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"  wrote {OUT}/{outstem}.png and .pdf")

def main():
    data = load()
    cols = [c for c in data["F0"].columns]
    means = [c for c in cols if c.endswith("_mean")
             and not any(m in c for m in ("vel", "acc", "divergence", "convergence",
                                          "vergere", "directionality", "node_"))]
    sums = [c for c in cols if c.endswith("_sum")
            and not any(m in c for m in ("vel", "acc", "divergence", "convergence",
                                         "vergere", "directionality", "node_"))]
    print(f"{len(means)} intensive, {len(sums)} extensive parameters")

    make_figure(
        data, means,
        "Mitochondrial network morphometry over 24 h of palmitate — shape parameters",
        "AdMSC, 250 µM palmitate added after the 0 h frame. Per-object means: these describe "
        "the organelles that were segmented, and are insensitive to how much of the network was "
        "captured.",
        "fig1_shape_parameters")

    make_figure(
        data, sums,
        "Mitochondrial network morphometry over 24 h of palmitate — summed parameters",
        "Totals over all objects in the field. These scale with how much network was captured "
        "and carry markedly more technical variance — read the trend, not the value. "
        "† = sum of a shape ratio, which reports object count more than size.",
        "fig2_summed_parameters", mark_ratio=True)

if __name__ == "__main__":
    main()
