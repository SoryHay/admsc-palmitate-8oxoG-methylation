#!/usr/bin/env bash
# Run Nellie on every time point of the 24 h time-lapse.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail

PROJ=/data/imaging_nellie
DATA=<DATA_DISK>/imaging_AdMSC_PA
PY="$PROJ/env/bin/python"
STEM="AdMSCs-24h-40x-250uMSP_2026-01-20_AdMSCs-24h-250uM_13.03.35"

export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=0

LINE="${1:-deconvolved}"
case "$LINE" in
  deconvolved) IMS_DIR="$DATA/raw/ims_deconvolved"; PREFIX="${STEM}_deconvolved" ;;
  raw)         IMS_DIR="$DATA/raw/ims_raw";         PREFIX="${STEM}" ;;
  *) echo "usage: $0 [deconvolved|raw]" >&2; exit 2 ;;
esac

echo "############ line: $LINE ############"

echo "===== [1] SUM-projection inputs (content arm) ====="
for F in F0 F1 F2; do
  "$PY" "$PROJ/scripts/ims_to_ometiff_mip.py" \
      "$IMS_DIR/${PREFIX}_${F}.ims" "$DATA/input_2Dt" --variant zsum
done

echo "===== [2] MAX-projection inputs, where missing (morphology arm) ====="
for F in F0 F1 F2; do
  if [ ! -f "$DATA/input_2Dt/${PREFIX}_${F}_mipwin1_2Dt.ome.tif" ]; then
    "$PY" "$PROJ/scripts/ims_to_ometiff_mip.py" \
        "$IMS_DIR/${PREFIX}_${F}.ims" "$DATA/input_2Dt" --variant mipwin
  else
    echo "  ${PREFIX}_${F}_mipwin1_2Dt.ome.tif already present"
  fi
done

echo "===== [3] Nellie, all 24 timepoints, one field at a time ====="
for F in F0 F1 F2; do
  OUT="$DATA/prod_${LINE}_${F}"
  echo "--- $F -> $OUT ---"
  rm -rf "$OUT"
  "$PY" "$PROJ/scripts/run_nellie.py" \
      "$DATA/input_2Dt/${PREFIX}_${F}_mipwin1_2Dt.ome.tif" "$OUT" \
      --device gpu --morphology-only
done

echo "===== DONE ($LINE) ====="
