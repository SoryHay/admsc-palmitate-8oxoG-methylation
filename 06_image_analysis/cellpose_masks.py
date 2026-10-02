"""Cellpose cell masks.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, os
import numpy as np, tifffile

ap = argparse.ArgumentParser()
ap.add_argument("--image", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--downscale", type=int, default=4)
ap.add_argument("--cellprob", type=float, default=-2.0)
ap.add_argument("--flow", type=float, default=0.8)
ap.add_argument("--model", default="cpsam_v2")
ap.add_argument("--stitch", type=float, default=0.0,
                help="stitch_threshold: segment the whole 2D+t stack at once and carry labels\n                     across time by overlap. Two independent trackers applied AFTER per-frame\n                     segmentation both returned only ~15 complete tracks out of ~50 cells,\n                     because per-frame instance boundaries are unstable - neighbouring cells\n                     merge in one frame and separate in the next. Stitching removes the linking\n                     step entirely instead of repairing it.")
a = ap.parse_args()

import torch
from cellpose import models
assert torch.cuda.is_available(), "GPU not visible in the container"
print(f"{a.model} on {torch.cuda.get_device_name(0)}, cellprob={a.cellprob} flow={a.flow}", flush=True)
m = models.CellposeModel(gpu=True, pretrained_model=a.model)

n_t = tifffile.TiffFile(a.image).series[0].shape[0]
os.makedirs(os.path.dirname(a.out), exist_ok=True)
out = []; frames = []
for t in range(n_t):
    fr = tifffile.imread(a.image, key=t)[::a.downscale, ::a.downscale].astype(np.float32)
    p1, p99 = np.percentile(fr, (1, 99.5))
    fr = np.clip((fr - p1) / max(p99 - p1, 1e-9), 0, 1)
    frames.append(fr)
if a.stitch > 0:
    stack = np.stack(frames)

    masks = m.eval(stack, cellprob_threshold=a.cellprob, flow_threshold=a.flow,
                   stitch_threshold=a.stitch, z_axis=0, batch_size=8)[0]
    masks = np.asarray(masks)
    out = [masks[t].astype(np.uint16) for t in range(masks.shape[0])]
    ids = [int(x.max()) for x in out]
    print(f"  stitched stack: {masks.max()} distinct labels overall, per-frame max {ids[:6]}...")
else:
    for t, fr in enumerate(frames):
        masks = m.eval(fr, cellprob_threshold=a.cellprob, flow_threshold=a.flow, batch_size=8)[0]
        out.append(masks.astype(np.uint16))
        print(f"  t{t:02d}: {int(masks.max()):3d} cells", flush=True)
np.save(a.out, np.stack(out))
print("saved", a.out)
