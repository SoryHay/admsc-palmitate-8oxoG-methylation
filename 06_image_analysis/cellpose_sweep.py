"""Choice of the Cellpose operating point.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import argparse, numpy as np, tifffile

ap = argparse.ArgumentParser()
ap.add_argument("--image", required=True)
ap.add_argument("--labels", required=True)
ap.add_argument("--t", type=int, default=0)
ap.add_argument("--downscale", type=int, default=4)
a = ap.parse_args()

import torch
from cellpose import models
print(f"torch {torch.__version__} cuda={torch.cuda.is_available()}", flush=True)

frame = tifffile.imread(a.image, key=a.t)
org = tifffile.imread(a.labels, key=a.t)
omask = org > 0
areas = np.bincount(org.ravel())[1:]
labs = np.unique(org)[1:]
XY = 0.15478489533094372

small = frame[::a.downscale, ::a.downscale].astype(np.float32)
p1, p99 = np.percentile(small, (1, 99.5))
small = np.clip((small - p1) / max(p99 - p1, 1e-9), 0, 1)

def score(masks):
    cells = np.repeat(np.repeat(masks, a.downscale, 0), a.downscale, 1)[:org.shape[0], :org.shape[1]]
    inside = omask & (cells > 0)
    nc = int(cells.max())
    o = org[omask].astype(np.int64); c = cells[omask].astype(np.int64)
    cnt = np.bincount(o * (nc + 1) + c, minlength=(int(org.max()) + 1) * (nc + 1))
    win = cnt.reshape(int(org.max()) + 1, nc + 1).argmax(axis=1)
    orph = win[labs] == 0
    oa = areas[labs[orph] - 1].mean() * XY**2 if orph.any() else float("nan")
    ka = areas[labs[~orph] - 1].mean() * XY**2
    return (inside.sum() / omask.sum() * 100, orph.mean() * 100, oa / ka,
            (cells > 0).sum() / cells.size * 100, nc)

print(f"\n{'model':12s} {'cellprob':>9s} {'flow':>5s} | {'inside%':>8s} {'orphan%':>8s} "
      f"{'bias':>6s} {'cover%':>7s} {'n':>4s}")
best = None
for model_name in ("cpsam_v2", "cpsam", "cpdino"):
    try:
        m = models.CellposeModel(gpu=True, pretrained_model=model_name)
    except Exception as e:
        print(f"{model_name:12s} unavailable: {str(e)[:60]}")
        continue
    for cp in (0.0, -2.0, -4.0):
        for fl in (0.4, 0.8):
            try:
                masks = m.eval(small, cellprob_threshold=cp, flow_threshold=fl, batch_size=8)[0]
            except Exception as e:
                print(f"{model_name:12s} {cp:9.1f} {fl:5.1f} | FAILED {str(e)[:40]}")
                continue
            s = score(masks)
            flag = ""
            if s[0] > 90 and s[3] < 55:
                flag = "  <--"
                if best is None or s[0] > best[1][0]:
                    best = ((model_name, cp, fl), s)
            print(f"{model_name:12s} {cp:9.1f} {fl:5.1f} | {s[0]:7.1f}% {s[1]:7.1f}% "
                  f"{s[2]:6.2f} {s[3]:6.1f}% {s[4]:4d}{flag}", flush=True)
print(f"\nbest by inside% subject to cover% < 55: {best}")
