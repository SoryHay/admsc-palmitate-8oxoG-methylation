#!/usr/bin/env bash
# Pooled 5mC per region for every region class and gene panel, both arms.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
P=/data/8oxo_project; B=$P/AdMSC_DMR/trackB_rerio_allctx; F=$P/AdMSC_DMR/features
R=$P/AdMSC_DMR/results; L=$P/AdMSC_DMR/logs; T=${CLAUDE_JOB_DIR:-/tmp}/tmp/l1
mkdir -p "$R" "$T"; st(){ date '+%H:%M:%S'; }

for arm in NT PA; do
  [ -s "$T/$arm.cpg.bed" ] || { echo "[$(st)] prep $arm"; \
    zcat "$B/$arm/ctx_thr070.cpg.tsv.gz" | tail -n +2 \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$6}' | sort -k1,1 -k2,2n > "$T/$arm.cpg.bed"; }
  echo "[$(st)] $arm CpG positions: $(wc -l < "$T/$arm.cpg.bed")"
done

for set in islands shores shelves promoters_tss1kb genebodies exons_merged panel_exons innate_exons first_exon first_intron promoters_by_cpg_class; do
  bed="$F/$set.bed"; [ -s "$bed" ] || continue
  echo "[$(st)] === $set ($(wc -l < "$bed") regions) ==="
  for arm in NT PA; do
    sort -k1,1 -k2,2n "$bed" | cut -f1-4 \
      | bedtools map -a - -b "$T/$arm.cpg.bed" -c 4,5 -o sum,sum -null 0 > "$T/$set.$arm.mapped.tsv"
  done
  python3 "$P/AdMSC_DMR/scripts/layer1_cpg_aggregate.py" \
     "$T/$set.NT.mapped.tsv" "$T/$set.PA.mapped.tsv" "$R/layer1_cpg_$set.tsv" 20 \
     2>&1 | tee "$L/layer1_$set.log"
done
echo "[$(st)] === LAYER 1 (own aggregation) COMPLETE ==="
