#!/usr/bin/env bash
# The same pooled aggregation at call thresholds 0.90 and 0.98.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
P=/data/8oxo_project; B=$P/AdMSC_DMR/trackB_rerio_allctx; F=$P/AdMSC_DMR/features
R=$P/AdMSC_DMR/results; L=$P/AdMSC_DMR/logs; T=${CLAUDE_JOB_DIR:-/tmp}/tmp/l1sw
mkdir -p "$R" "$T"; st(){ date '+%H:%M:%S'; }
for thr in 090 098; do
  for arm in NT PA; do
    o="$B/$arm/ctxsw_thr$thr"
    [ -s "$o.cpg.tsv.gz" ] || { echo "[$(st)] classify $arm thr$thr"; \
      python3 "$P/AdMSC_DMR/scripts/layer_classify_context.py" \
        "$B/$arm/pileup_feat_thr$thr.bed.gz" "$F/union_features.bed" \
        "$P/refs/T2T-CHM13v2.0.fa" "$o" 2>&1 | tail -3; }
    [ -s "$T/$arm.$thr.bed" ] || zcat "$o.cpg.tsv.gz" | tail -n +2 \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$6}' | sort -k1,1 -k2,2n > "$T/$arm.$thr.bed"
  done
  for set in islands promoters_tss1kb genebodies shores; do
    for arm in NT PA; do
      sort -k1,1 -k2,2n "$F/$set.bed" | cut -f1-4 \
        | bedtools map -a - -b "$T/$arm.$thr.bed" -c 4,5 -o sum,sum -null 0 > "$T/$set.$arm.$thr.tsv"
    done
    echo "[$(st)] --- $set @ 0.$thr ---"
    python3 "$P/AdMSC_DMR/scripts/layer1_cpg_aggregate.py" \
      "$T/$set.NT.$thr.tsv" "$T/$set.PA.$thr.tsv" "$R/layer1_sweep_${set}_thr$thr.tsv" 20 \
      2>&1 | tee "$L/layer1_sweep_${set}_$thr.log" | grep -E "pooled|regions tested"
  done
done
echo "[$(st)] === SWEEP COMPLETE ==="
