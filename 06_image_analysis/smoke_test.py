"""Environment check.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os
import sys

print("=" * 62)
print("CUDA_VISIBLE_DEVICES =", os.environ.get("CUDA_VISIBLE_DEVICES", "<unset>"))
print("=" * 62)

import cupy as cp

n_dev = cp.cuda.runtime.getDeviceCount()
print(f"[cupy] version           : {cp.__version__}")
print(f"[cupy] runtime CUDA      : {cp.cuda.runtime.runtimeGetVersion()}")
print(f"[cupy] driver CUDA       : {cp.cuda.runtime.driverGetVersion()}")
print(f"[cupy] visible devices   : {n_dev}")
if n_dev != 1:
    print(f"  !! expected exactly 1 visible device, got {n_dev} — the pin did not take")

props = cp.cuda.runtime.getDeviceProperties(0)
free_b, total_b = cp.cuda.runtime.memGetInfo()
print(f"[cupy] device 0          : {props['name'].decode()}")
print(f"[cupy] VRAM free/total   : {free_b / 2**30:.1f} / {total_b / 2**30:.1f} GiB")

import cupyx.scipy.ndimage as cndi

a = cp.random.randint(0, 4096, size=(6, 2040, 2040), dtype=cp.uint16).astype(cp.float32)
g = cndi.gaussian_filter(a, sigma=(0, 2, 2))
m = cndi.maximum_filter(g, size=(1, 3, 3))
lbl, n_obj = cndi.label(g > g.mean() + 2 * g.std())
cp.cuda.Stream.null.synchronize()

peak = cp.get_default_memory_pool().used_bytes()
print(f"[gpu ] gaussian+maximum+label on a real-sized frame: OK, {int(n_obj)} objects")
print(f"[gpu ] pool in use       : {peak / 2**20:.0f} MiB "
      f"({peak / total_b * 100:.1f} % of card)")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import napari
print(f"[napari] version         : {napari.__version__}")

from npe2 import PluginManager

pm = PluginManager.instance()
pm.discover()
names = sorted(m.name for m in pm.iter_manifests())
print(f"[npe2] plugins found     : {names}")
print(f"[npe2] nellie registered : {'nellie' in names}")

import nellie
from nellie.run import run
from nellie.tracking.flow_interpolation import FlowInterpolator

print(f"[nellie] version         : {getattr(nellie, '__version__', 'n/a')}")
print(f"[nellie] device_type seen : {nellie.device_type}")
print(f"[nellie] run/FlowInterpolator importable: {callable(run) and FlowInterpolator is not None}")

print("=" * 62)
print("SMOKE TEST PASSED")
sys.exit(0)
