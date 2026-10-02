#!/usr/bin/env bash
# Install the image-analysis environment.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail

PROJ=/data/imaging_nellie
ENV_PREFIX="$PROJ/env"
PY_VERSION=3.11
CUPY_PKG=cupy-cuda12x
EXTRA_PKGS="h5py tifffile ome-types"

mkdir -p "$PROJ/logs" "$PROJ/scripts"

echo "== [1/3] conda env: python $PY_VERSION at $ENV_PREFIX =="
conda create -y -p "$ENV_PREFIX" -c conda-forge "python=$PY_VERSION"

PIP="$ENV_PREFIX/bin/pip"
echo "== [2/3] pip: napari (latest) + nellie + $EXTRA_PKGS =="
"$PIP" install --no-input "napari[all]" nellie $EXTRA_PKGS

echo "== [3/3] pip: $CUPY_PKG (GPU acceleration) =="
"$PIP" install --no-input "$CUPY_PKG"

echo "== versions =="
"$ENV_PREFIX/bin/python" - <<'PY'
import importlib.metadata as md
for pkg in ("napari", "nellie", "cupy-cuda12x", "numpy", "scikit-image", "nd2", "PyQt6"):
    try:
        print(f"{pkg:16s} {md.version(pkg)}")
    except md.PackageNotFoundError:
        print(f"{pkg:16s} NOT INSTALLED")
PY

echo "DONE"
