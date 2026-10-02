"""Imaris .ims to OME-TIFF with pixel calibration; maximum-intensity projections.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime

import h5py
import numpy as np
import tifffile

CAMERA_PIXEL_UM = 6.5
MAGNIFICATION = 40.0
MAG_CORRECTION = 1.049844041
XY_UM = CAMERA_PIXEL_UM / MAGNIFICATION / MAG_CORRECTION

def _s(v):
    """Imaris stores attributes as char arrays."""
    if isinstance(v, np.ndarray):
        return "".join(c.decode() for c in v)
    return v

def read_geometry(f):
    """Return (X, Y, Z, z_um, timestamps) from the Imaris metadata."""
    a = {k: _s(v) for k, v in f["DataSetInfo/Image"].attrs.items()}
    X, Y, Z = int(a["X"]), int(a["Y"]), int(a["Z"])
    z_um = (float(a["ExtMax2"]) - float(a["ExtMin2"])) / Z

    ta = {k: _s(v) for k, v in f["DataSetInfo/TimeInfo"].attrs.items()}
    keys = sorted((k for k in ta if k.startswith("TimePoint")),
                  key=lambda x: int(x.replace("TimePoint", "")))
    stamps = []
    for k in keys:
        try:
            stamps.append(datetime.strptime(ta[k].strip(), "%Y-%m-%d %H:%M:%S.%f"))
        except ValueError:
            stamps.append(None)
    return X, Y, Z, z_um, stamps

def tenengrad(plane: np.ndarray) -> float:
    """Standard autofocus measure: mean squared gradient magnitude. Higher = better focused."""
    gy, gx = np.gradient(plane.astype(np.float32))
    return float(np.mean(gx * gx + gy * gy))

def tile_focus(plane: np.ndarray, n: int = 3) -> np.ndarray:
    """Tenengrad on an n x n grid, to detect sample tilt (focus varying across the field)."""
    h, w = plane.shape
    out = np.empty((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(n):
            out[i, j] = tenengrad(plane[i * h // n:(i + 1) * h // n,
                                        j * w // n:(j + 1) * w // n])
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ims")
    ap.add_argument("out_dir")
    ap.add_argument("--window", type=int, default=1,
                    help="half-width in planes for the mipwin variant (default 1 -> 3 planes)")
    ap.add_argument("--variant", choices=["focus", "mipwin", "zsum", "all"], default="all",
                    help="focus = best plane only; mipwin = MAX over the drift-following window "
                         "(morphology); zsum = SUM over the same window (content); all = every one")
    ap.add_argument("--tiles", type=int, default=3, help="grid for the tilt check")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.ims))[0]

    with h5py.File(args.ims, "r") as f:
        X, Y, Z, z_um, stamps = read_geometry(f)
        res = f["DataSet/ResolutionLevel 0"]
        tps = sorted((k for k in res if k.startswith("TimePoint")),
                     key=lambda x: int(x.split()[-1]))
        n_t = len(tps)

        print(f"{base}: {n_t} timepoints, {Z} planes, {Y}x{X} px")
        print(f"  voxel XY : {XY_UM:.6f} um  (file said {(float(_s(f['DataSetInfo/Image'].attrs['ExtMax0'])) - float(_s(f['DataSetInfo/Image'].attrs['ExtMin0']))) / X:.6f}; "
              f"correction factor {CAMERA_PIXEL_UM / (CAMERA_PIXEL_UM / 6):.0f}x)")
        print(f"  voxel Z  : {z_um:.4f} um   (kept as-is; the Z extent matches the protocol)")

        deltas = [(stamps[i + 1] - stamps[i]).total_seconds()
                  for i in range(len(stamps) - 1) if stamps[i] and stamps[i + 1]]
        dt = float(np.median(deltas)) if deltas else float("nan")
        print(f"  interval : median {dt:.2f} s over {len(deltas)} gaps "
              f"(min {min(deltas):.2f}, max {max(deltas):.2f})")

        focus_rows, best_planes, tilt_rows = [], [], []
        stack_focus = np.empty((n_t, Y, X), dtype=np.uint16)
        stack_mip = np.empty((n_t, Y, X), dtype=np.uint16)

        stack_sum = np.empty((n_t, Y, X), dtype=np.uint32)

        for t, key in enumerate(tps):
            vol = res[f"{key}/Channel 0/Data"][:Z, :Y, :X]
            scores = np.array([tenengrad(vol[z]) for z in range(Z)])
            best = int(np.argmax(scores))
            best_planes.append(best)
            focus_rows.append([t] + [f"{s:.4f}" for s in scores] + [str(best)])

            tile_scores = np.stack([tile_focus(vol[z], args.tiles) for z in range(Z)])
            per_tile_best = tile_scores.argmax(axis=0)
            tilt_rows.append([t, best, per_tile_best.min(), per_tile_best.max(),
                              per_tile_best.max() - per_tile_best.min()])

            stack_focus[t] = vol[best]
            lo, hi = max(0, best - args.window), min(Z, best + args.window + 1)
            stack_mip[t] = vol[lo:hi].max(axis=0)
            stack_sum[t] = vol[lo:hi].astype(np.uint32).sum(axis=0)
            print(f"  T{t:02d}: best plane Z{best}  tile-best range Z{per_tile_best.min()}-Z{per_tile_best.max()}  "
                  f"mip window Z{lo}-Z{hi - 1}")

    meta = {
        "axes": "TYX",
        "PhysicalSizeX": XY_UM, "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": XY_UM, "PhysicalSizeYUnit": "µm",
        "TimeIncrement": dt, "TimeIncrementUnit": "s",
    }
    written = []
    if args.variant in ("focus", "all"):
        p = os.path.join(args.out_dir, f"{base}_focusplane_2Dt.ome.tif")
        tifffile.imwrite(p, stack_focus, ome=True, photometric="minisblack", metadata=meta)
        written.append(p)
    if args.variant in ("mipwin", "all"):
        p = os.path.join(args.out_dir, f"{base}_mipwin{args.window}_2Dt.ome.tif")
        tifffile.imwrite(p, stack_mip, ome=True, photometric="minisblack", metadata=meta)
        written.append(p)
    if args.variant in ("zsum", "all"):
        p = os.path.join(args.out_dir, f"{base}_zsum{args.window}_2Dt.ome.tif")
        tifffile.imwrite(p, stack_sum, ome=True, photometric="minisblack", metadata=meta)
        written.append(p)

    with open(os.path.join(args.out_dir, f"{base}_focus_profile.tsv"), "w") as fh:
        fh.write("t\t" + "\t".join(f"tenengrad_Z{z}" for z in range(Z)) + "\tbest_plane\n")
        for r in focus_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(args.out_dir, f"{base}_tilt_check.tsv"), "w") as fh:
        fh.write("t\tbest_plane_whole_field\ttile_best_min\ttile_best_max\ttile_spread\n")
        for r in tilt_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(args.out_dir, f"{base}_conversion.json"), "w") as fh:
        json.dump({
            "source_ims": os.path.abspath(args.ims),
            "outputs": written,
            "n_timepoints": n_t, "n_planes_in_source": Z, "shape_yx": [Y, X],
            "voxel_xy_um_corrected": XY_UM,
            "voxel_xy_um_as_recorded": XY_UM / 6.0,
            "correction_factor": 6.0,
            "correction_derivation": ("camera pixel 6.5 um = 12240 sub-units x 1.0833333 um / 2040 px; "
                                      "Imaris divided the SUB-UNIT by mag 40 x 1.049844 instead of the pixel"),
            "voxel_z_um": z_um,
            "time_interval_s_median": dt,
            "time_interval_s_min": min(deltas) if deltas else None,
            "time_interval_s_max": max(deltas) if deltas else None,
            "timestamps": [s.isoformat() if s else None for s in stamps],
            "best_focus_plane_per_t": best_planes,
            "mip_window_half_width": args.window,
        }, fh, indent=2)

    print("\nwrote:")
    for p in written:
        print("  ", p, f"{os.path.getsize(p) / 2**20:.0f} MiB")

if __name__ == "__main__":
    main()
