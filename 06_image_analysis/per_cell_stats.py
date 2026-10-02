"""Attribute mitochondrial objects to cells; per-cell statistics (Figure S3).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations

import argparse, glob, os
import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from skimage.morphology import disk

def _within_cell(a):
    """
    Statistics of the organelle-size distribution INSIDE one cell.

    A per-cell mean cannot distinguish a uniformly shrunken network from one in which a large
    connected structure broke up while the rest was untouched - and fission predicts the second.
    These describe the shape of the distribution instead of collapsing it to its mean.
    """
    a = np.asarray(a, dtype=float)
    a = a[a > 0]
    if a.size == 0:
        return {}
    tot = a.sum()
    return {
        "median_organelle_area_um2": float(np.median(a)),
        "p90_organelle_area_um2": float(np.percentile(a, 90)),

        "area_weighted_size_um2": float((a ** 2).sum() / tot),
        "frac_area_largest": float(a.max() / tot),
        "cv_organelle_area": float(a.std(ddof=1) / a.mean()) if a.size > 1 else float("nan"),
    }

ap = argparse.ArgumentParser()
ap.add_argument("--field", required=True)
ap.add_argument("--prod-dir", required=True)
ap.add_argument("--masks", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--downscale", type=int, default=4)
ap.add_argument("--xy-um", type=float, default=0.15478489533094372)
ap.add_argument("--closing-um", type=float, default=4.0)
a = ap.parse_args()
XY, PX2 = a.xy_um, a.xy_um ** 2

inst = glob.glob(os.path.join(a.prod_dir, "nellie_necessities", "*im_instance_label.ome.tif"))[0]
feat = pd.read_csv(glob.glob(os.path.join(a.prod_dir, "*features_organelles.csv"))[0])

_DROP = ("vel", "acc", "divergence", "convergence", "vergere", "directionality", "node_")
_COORD = ("x_raw", "y_raw", "z_raw", "reassigned_label_raw")
keep = [c for c in feat.columns if c in ("t", "label") or
        (not any(m in c for m in _DROP) and c not in _COORD
         and c.endswith(("_mean", "_sum", "_raw")))]
feat = feat[keep]
cp = np.load(a.masks)
n_t = tifffile.TiffFile(inst).series[0].shape[0]
r = max(1, int(round(a.closing_um / XY)))

rows, maprows = [], []
for t in range(n_t):
    org = tifffile.imread(inst, key=t)
    cells = np.repeat(np.repeat(cp[t], a.downscale, 0), a.downscale, 1)[:org.shape[0], :org.shape[1]]
    omask = org > 0

    territory = ndi.binary_fill_holes(ndi.binary_closing(omask, disk(r)))

    n_cells = int(cells.max())
    o = org[omask].astype(np.int64); c = cells[omask].astype(np.int64)
    cnt = np.bincount(o * (n_cells + 1) + c,
                      minlength=(int(org.max()) + 1) * (n_cells + 1)).reshape(int(org.max()) + 1, n_cells + 1)
    owner = cnt.argmax(axis=1)
    labs = np.unique(org)[1:]
    for l in labs:
        maprows.append((t, int(l), int(owner[l])))

    border = set(np.unique(np.concatenate([cells[0], cells[-1], cells[:, 0], cells[:, -1]]))) - {0}
    org_area = np.bincount(org.ravel(), minlength=int(org.max()) + 1) * PX2
    cell_px = np.bincount(cells.ravel(), minlength=n_cells + 1)
    terr_px = np.bincount(cells[territory].ravel(), minlength=n_cells + 1)

    own = owner[labs]
    for cid in range(1, n_cells + 1):
        mine = labs[own == cid]
        if len(mine) == 0:
            continue
        area_um2 = terr_px[cid] * PX2
        rows.append({
            "field": a.field, "t": t, "cell_id": cid,
            "touches_border": cid in border,
            "cell_area_um2": area_um2,
            "cellpose_area_um2": cell_px[cid] * PX2,
            "n_organelles": len(mine),
            "mito_area_um2": org_area[mine].sum(),
            "mito_area_fraction": org_area[mine].sum() / area_um2 if area_um2 > 0 else np.nan,
            "organelle_density_per_100um2": len(mine) / area_um2 * 100 if area_um2 > 0 else np.nan,
            "mean_organelle_area_um2": org_area[mine].mean(),
            **_within_cell(org_area[mine]),
        })
    orph = (own == 0).sum()
    print(f"  {a.field} t{t:02d}: {n_cells:3d} cells, {len(labs):5d} organelles, "
          f"{orph:4d} orphan ({orph/max(len(labs),1)*100:4.1f} %)", flush=True)

per_cell = pd.DataFrame(rows)
mp = pd.DataFrame(maprows, columns=["t", "label", "cell_id"])
merged = feat.merge(mp, on=["t", "label"], how="left")
shape_cols = [c for c in feat.columns if c.endswith(("_mean", "_raw"))]
if shape_cols:
    agg = merged[merged.cell_id > 0].groupby(["t", "cell_id"])[shape_cols].mean().reset_index()
    per_cell = per_cell.merge(agg, on=["t", "cell_id"], how="left")

os.makedirs(os.path.dirname(a.out), exist_ok=True)
per_cell.to_csv(a.out, index=False)
print(f"\nwrote {len(per_cell)} cell-timepoint rows -> {a.out}")
