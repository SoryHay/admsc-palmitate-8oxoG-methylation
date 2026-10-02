"""Nellie mitochondrial network segmentation at a fixed operating point.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
from __future__ import annotations

import argparse
import json
import os
import time

MAX_MOVE_UM = 10.0
REMOVE_EDGES = False

OTSU_THRESH = False
THRESHOLD = None

def patch_flow_interpolator(max_move_um: float):
    """Cap the voxel-matching radius at `max_move_um` microns regardless of dt."""
    import nellie.tracking.flow_interpolation as fi_mod
    import nellie.tracking.voxel_reassignment as vr_mod
    import nellie.feature_extraction.hierarchical as hi_mod

    base = fi_mod.FlowInterpolator

    class GuardedFlowInterpolator(base):
        def __init__(self, im_info, num_t=None, max_distance_um=None, forward=True):
            if max_distance_um is None:
                dt = float(im_info.dim_res["T"])
                max_distance_um = max_move_um / dt if dt > 0 else 0.5
            super().__init__(im_info, num_t=num_t,
                             max_distance_um=max_distance_um, forward=forward)

    for mod in (fi_mod, vr_mod, hi_mod):
        mod.FlowInterpolator = GuardedFlowInterpolator
    return GuardedFlowInterpolator

def patch_hu_tracking_host_coercion():
    """
    Work around an upstream bug in nellie 1.0.4.

    `HuMomentTracking._match_frames_sparse` (hu_tracking.py:947) declares in its own docstring that
    it "operates entirely on CPU/NumPy to reduce GPU memory pressure", and calls np.asarray() on its
    arguments — but under the GPU backend it is handed cupy arrays for `stats_vecs` / `hu_vecs`,
    and cupy refuses implicit conversion:

        TypeError: Implicit conversion to a NumPy array is not allowed. Please use `.get()`

    The sparse branch is taken whenever the frame has many objects, which is always the case here
    (~50k objects per frame), so the GPU path cannot complete without this.

    The patch only moves device arrays to the host, which is what the function already intends.
    It changes no numerical behaviour. Reported here so it is not mistaken for a tuning choice.
    """
    import numpy as np
    import nellie.tracking.hu_tracking as hu_mod

    def to_host(x):
        if hasattr(x, "get") and type(x).__module__.startswith("cupy"):
            return x.get()
        if isinstance(x, (list, tuple)):
            return type(x)(to_host(v) for v in x)
        return x

    original = hu_mod.HuMomentTracking._match_frames_sparse

    def patched(self, *args, **kwargs):
        return original(self, *(to_host(a) for a in args),
                        **{k: to_host(v) for k, v in kwargs.items()})

    hu_mod.HuMomentTracking._match_frames_sparse = patched
    _ = np

def run_morphology_only(file_info, device, timeit=True):
    """
    Filter -> Label -> Network -> Markers -> Hierarchy(no motility, no nodes).

    `nellie.run.run()` always runs six stages and constructs Hierarchy with `skip_nodes=False`,
    overriding the class's own `skip_nodes=True` default. On this dataset that is expensive for
    nothing: of 2365 s for four timepoints, segmentation took 5.5 s, HuMomentTracking 466 s,
    VoxelReassigner 387 s and Hierarchy 1507 s — i.e. 99.8 % of the time went into the tracking
    and node families, every one of which we have already established is not reportable at a
    3600 s frame interval.

    Skipping is safe, not a hack:
      - Hierarchy loads the reassigned-label memmaps only `if os.path.exists(...)`
        (hierarchical.py:224-228), so it degrades cleanly when VoxelReassigner never ran;
      - with `enable_motility=False` it never builds the flow interpolators
        (hierarchical.py:541-552) and fills the motion columns with NaN by design.

    The morphology and topology columns are computed by the same code either way — which is
    checked, not assumed, by comparing this path against the full run on the same input.
    """
    import time as _time

    from nellie.feature_extraction.hierarchical import Hierarchy
    from nellie.im_info.verifier import ImInfo
    from nellie.segmentation.filtering import Filter
    from nellie.segmentation.labelling import Label
    from nellie.segmentation.mocap_marking import Markers
    from nellie.segmentation.networking import Network

    im_info = ImInfo(file_info)
    timings = {}
    for name, stage in (
        ("Filter", lambda: Filter(im_info, remove_edges=REMOVE_EDGES, device=device)),
        ("Label", lambda: Label(im_info, otsu_thresh_intensity=OTSU_THRESH,
                                threshold=THRESHOLD, device=device)),
        ("Network", lambda: Network(im_info, device=device)),
        ("Markers", lambda: Markers(im_info, device=device)),
        ("Hierarchy", lambda: Hierarchy(im_info, skip_nodes=True, enable_motility=False,
                                        enable_adjacency=False, device=device)),
    ):
        t0 = _time.perf_counter()
        stage().run()
        timings[name] = _time.perf_counter() - t0
        if timeit:
            print(f"Nellie (morphology-only): {name} step took {timings[name]:.4f} seconds")
    return im_info, timings

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output_dir")
    ap.add_argument("--t-start", type=int, default=0)
    ap.add_argument("--t-end", type=int, default=None)
    ap.add_argument("--device", default="gpu", choices=["auto", "cpu", "gpu"])
    ap.add_argument("--morphology-only", action="store_true",
                    help="skip the tracking stages and the node/motility feature families "
                         "(they are not reportable at dt=3600 s); see run_morphology_only()")
    args = ap.parse_args()

    from nellie.im_info.verifier import FileInfo
    from nellie.run import run

    patch_flow_interpolator(MAX_MOVE_UM)
    patch_hu_tracking_host_coercion()

    fi = FileInfo(args.input, output_dir=args.output_dir)
    fi.find_metadata()
    fi.load_metadata()
    print(f"axes   : {fi.axes}")
    print(f"shape  : {fi.shape}")
    print(f"dim_res: {fi.dim_res}")

    if args.t_end is not None or args.t_start:
        fi.select_temporal_range(args.t_start, args.t_end)
    fi._validate()
    print(f"good_axes={fi.good_axes} good_dims={fi.good_dims} "
          f"errors={fi.get_validation_errors()}")
    if not (fi.good_axes and fi.good_dims):
        raise SystemExit("metadata did not validate — refusing to run on unvalidated dimensions")

    dt = float(fi.dim_res["T"])
    print(f"guard rail: effective matching radius = {MAX_MOVE_UM} um "
          f"(nellie default would be {0.5 * dt:.0f} um at dt={dt:.1f} s)")

    t0 = time.perf_counter()
    if args.morphology_only:
        im_info, _ = run_morphology_only(fi, device=args.device)
    else:
        im_info = run(fi, remove_edges=REMOVE_EDGES, otsu_thresh_intensity=OTSU_THRESH,
                      threshold=THRESHOLD, timeit=True, device=args.device)
    elapsed = time.perf_counter() - t0

    print(f"\ntotal wall time: {elapsed:.1f} s")

    manifest = {
        "input": os.path.abspath(args.input),
        "output_dir": os.path.abspath(args.output_dir),
        "t_start": fi.t_start, "t_end": fi.t_end,
        "axes": fi.axes, "shape": list(fi.shape), "dim_res": fi.dim_res,
        "device": args.device,
        "morphology_only": args.morphology_only,
        "operating_point": {
            "max_move_um_guard_rail": MAX_MOVE_UM,
            "remove_edges": REMOVE_EDGES,
            "otsu_thresh_intensity": OTSU_THRESH,
            "threshold": THRESHOLD,
        },
        "wall_time_s": elapsed,
        "tracking_outputs_reportable": False,
        "tracking_note": ("dt=3600 s; no frame-to-frame organelle correspondence exists at this "
                          "sampling. Velocity/acceleration and voxel-reassignment outputs are "
                          "exploratory only."),
    }
    os.makedirs(args.output_dir, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.input))[0]
    with open(os.path.join(args.output_dir, f"{name}_run_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print("manifest written")

if __name__ == "__main__":
    main()
