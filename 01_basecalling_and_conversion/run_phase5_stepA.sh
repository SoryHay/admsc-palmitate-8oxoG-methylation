#!/usr/bin/env bash
# Split each pileup into CpG and CpH from the reference sequence.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
P=/data/8oxo_project; B=$P/AdMSC_DMR/trackB_rerio_allctx; F=$P/AdMSC_DMR/features
L=$P/AdMSC_DMR/logs; R=$P/refs/T2T-CHM13v2.0.fa; mkdir -p "$L" "$P/AdMSC_DMR/results"
st(){ date '+%H:%M:%S'; }
for arm in NT PA; do
  for thr in 070 095; do
    o="$B/$arm/ctx_thr$thr"
    [ -s "$o.cpg.tsv.gz" ] && { echo "[$(st)] $arm thr$thr present, skip"; continue; }
    echo "[$(st)] --- classify $arm thr$thr ---"
    python3 "$P/AdMSC_DMR/scripts/layer_classify_context.py" \
      "$B/$arm/pileup_gw_thr$thr.bed.gz" "$F/union_features.bed" "$R" "$o" \
      2>&1 | tee "$L/classify_${arm}_thr$thr.log"
  done
done
echo "[$(st)] === STEP A COMPLETE ==="
