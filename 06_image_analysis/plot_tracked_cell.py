"""Figure 6a.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations
import argparse, glob
import numpy as np, tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.colors as mc, matplotlib.gridspec as gridspec
from scipy import ndimage as ndi

ap = argparse.ArgumentParser()
ap.add_argument("--field", default="F0")
ap.add_argument("--cell", type=int, default=30)
ap.add_argument("--times", type=int, nargs="+", default=[0, 4, 8, 12, 16, 20])
ap.add_argument("--window-um", type=float, default=110.0)
ap.add_argument("--dpi", type=int, default=300)
ap.add_argument("--xy-um", type=float, default=0.15478489533094372)
ap.add_argument("--downscale", type=int, default=4)
ap.add_argument("--out", default="/data/imaging_nellie/figures/fig5_tracked_cell_network")
a = ap.parse_args()
XY, DS, F, L = a.xy_um, a.downscale, a.field, a.cell

D = "<DATA_DISK>/imaging_AdMSC_PA"
stem = "AdMSCs-24h-40x-250uMSP_2026-01-20_AdMSCs-24h-250uM_13.03.35_deconvolved"
img = f"{D}/input_2Dt/{stem}_{F}_mipwin1_2Dt.ome.tif"
inst = glob.glob(f"{D}/prod_deconvolved_{F}/nellie_necessities/*im_instance_label.ome.tif")[0]
skelp = glob.glob(f"{D}/prod_deconvolved_{F}/nellie_necessities/*im_skel.ome.tif")[0]
m = np.load(f"/data/imaging_nellie/masks/cellpose_stitched_{F}.npy")
n_t = m.shape[0]

cent = []
for t in range(n_t):
    sel = np.repeat(np.repeat(m[t] == L, DS, 0), DS, 1)
    cy, cx = ndi.center_of_mass(sel)
    cent.append((cx, cy))
cent = np.array(cent)

half = int(a.window_um / XY / 2); S = 2 * half
COL = np.array(mc.hsv_to_rgb((0.13, 0.9, 1.0)))
SKEL = np.array([0.15, 0.95, 0.85])

def crop(arr, cx, cy):
    out = np.zeros((S, S), arr.dtype)
    y0, y1, x0, x1 = int(cy) - half, int(cy) + half, int(cx) - half, int(cx) + half
    sub = arr[max(0, y0):min(arr.shape[0], y1), max(0, x0):min(arr.shape[1], x1)]
    out[max(0, -y0):max(0, -y0) + sub.shape[0], max(0, -x0):max(0, -x0) + sub.shape[1]] = sub
    return out

fig = plt.figure(figsize=(20.5, 4.3))
grid = gridspec.GridSpec(1, len(a.times) + 2,
                         width_ratios=[1] * len(a.times) + [0.12, 1.55], wspace=0.04)
for i, t in enumerate(a.times):
    ax = fig.add_subplot(grid[0, i]); cx, cy = cent[t]
    base = crop(tifffile.imread(img, key=t).astype(np.float32), cx, cy)
    lo, hi = np.percentile(base[base > 0], (2, 99.8)) if (base > 0).any() else (0, 1)
    base = np.clip((base - lo) / (hi - lo + 1e-9), 0, 1)
    cell = crop(np.repeat(np.repeat(m[t] == L, DS, 0), DS, 1).astype(np.uint8), cx, cy).astype(bool)
    org = crop((tifffile.imread(inst, key=t) > 0).astype(np.uint8), cx, cy).astype(bool) & cell
    sk = crop((tifffile.imread(skelp, key=t) > 0).astype(np.uint8), cx, cy).astype(bool) & cell
    rgb = np.stack([base] * 3, -1) * 0.75
    rgb[org] = np.clip(rgb[org] * 0.30 + COL * 0.70, 0, 1)
    rgb[ndi.binary_dilation(sk, np.ones((2, 2)))] = SKEL
    rgb[cell ^ ndi.binary_erosion(cell, np.ones((7, 7)))] = [1, 1, 1]
    ax.imshow(rgb, interpolation="nearest"); ax.axis("off")
    ax.set_title(f"{t} h", fontsize=13, color="#22292d")
    if i == 0:
        n = 25 / XY
        ax.plot([14, 14 + n], [S - 18] * 2, color="w", lw=3)
        ax.text(14 + n / 2, S - 26, "25 µm", color="w", ha="center", va="bottom", fontsize=9)

ax = fig.add_subplot(grid[0, -1])
base = tifffile.imread(img, key=0)[::DS, ::DS].astype(np.float32)
lo, hi = np.percentile(base, (2, 99.7)); base = np.clip((base - lo) / (hi - lo), 0, 1)
rgb = np.stack([base] * 3, -1) * 0.8
c0 = m[0] == L
rgb[c0 ^ ndi.binary_erosion(c0, np.ones((3, 3)))] = COL
ax.imshow(rgb, interpolation="nearest"); ax.axis("off")
tr = cent / DS
sc = ax.scatter(tr[:, 0], tr[:, 1], c=np.arange(n_t), cmap="autumn_r", s=16, zorder=4)
ax.plot(tr[:, 0], tr[:, 1], color="#ff3b1f", lw=1.4, alpha=0.9, zorder=3)
net = np.hypot(*(cent[-1] - cent[0])) * XY
path = np.sum(np.hypot(*np.diff(cent, axis=0).T)) * XY
ax.set_title(f"migration over 23 h — net {net:.0f} µm, path {path:.0f} µm", fontsize=11, color="#22292d")
cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
cb.set_label("hours", fontsize=8); cb.ax.tick_params(labelsize=7)

fig.suptitle("A single tracked AdMSC through 24 h of palmitate — white, cell outline; yellow, its "
             "segmented mitochondria; cyan, the skeletonised network; right, its migration track",
             fontsize=12.5, fontweight="bold", color="#22292d", y=1.03)
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=a.dpi, bbox_inches="tight", facecolor="white")
print(f"wrote {a.out}.png at {a.dpi} dpi  (net {net:.0f} um, path {path:.0f} um)")
