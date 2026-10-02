#!/usr/bin/env bash
# Alignment and modkit pileup at the declared call thresholds.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
PROJ=/data/8oxo_project
BASE=$PROJ/AdMSC_DMR/trackB_rerio_allctx
FEAT=$PROJ/AdMSC_DMR/features
LOG=$PROJ/AdMSC_DMR/logs
REF=$PROJ/refs/T2T-CHM13v2.0.fa
MODKIT=$(ls /data/modkit/latest/*/modkit | head -1)
THREADS=12
cd "$PROJ"
stamp(){ date '+%Y-%m-%d %H:%M:%S'; }

if [ ! -s "$FEAT/union_features.bed" ]; then
  cat "$FEAT"/{islands,shores,shelves,promoters_tss1kb,genebodies,exons_merged,panel_exons,innate_exons}.bed 2>/dev/null \
    | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > "$FEAT/union_features.bed"
fi
echo "[$(stamp)] union feature BED: $(wc -l < "$FEAT/union_features.bed") intervals, $(awk '{s+=$3-$2} END{printf "%.2f", s/1e9}' "$FEAT/union_features.bed") Gb"

for arm in NT PA; do
  D=$BASE/$arm
  if [ -s "$D/aln_T2T.bam.bai" ]; then echo "[$(stamp)] $arm alignment already present, skipping"; continue; fi
  echo "[$(stamp)] === PHASE 3 $arm — align $(ls "$D"/pass/*.bam | wc -l) pass BAMs to T2T ==="
  t0=$SECONDS
  bash "$PROJ/scripts/align_mods_t2t.sh" "$D" "$THREADS" 4 > "$LOG/align_$arm.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$D/aln_T2T.bam" ]; then
    echo "[$(stamp)] ⚠️ $arm ALIGNMENT FAILED (rc=$rc) — stopping"; exit 1
  fi
  PRIM=$(samtools view -c -F 0x904 -@ 4 "$D/aln_T2T.bam")
  TAGS=$(samtools view -F 0x904 "$D/aln_T2T.bam" 2>/dev/null | head -2000 | grep -c 'MM:Z:')
  echo "[$(stamp)] === $arm ALIGNED — $((( SECONDS-t0 )/60)) min, $PRIM primary, $TAGS/2000 carry MM:Z ==="
  if [ "$TAGS" -lt 1900 ]; then echo "[$(stamp)] ⚠️ $arm TAGS LOST — update-tags or -y failed. STOPPING"; exit 1; fi
  rm -f "$D/merged.bam" "$D/upd.bam"
  echo "[$(stamp)] $arm intermediates deleted (merged.bam, upd.bam)"
done

pileup(){    # $1=arm  $2=threshold  $3=tag  $4=optional --include-bed path
  local arm=$1 thr=$2 tag=$3 bed=${4:-}
  local out="$BASE/$arm/pileup_${tag}.bed.gz"
  [ -s "$out.tbi" ] && { echo "[$(stamp)] $arm $tag already present, skipping"; return 0; }
  local extra=(); [ -n "$bed" ] && extra=(--include-bed "$bed" --modified-bases C:m)
  echo "[$(stamp)] --- pileup $arm thr=$thr ${bed:+(features only)} ---"
  local t0=$SECONDS
  "$MODKIT" pileup "$BASE/$arm/aln_T2T.bam" - --ref "$REF" \
      --filter-threshold C:"$thr" --filter-threshold A:"$thr" \
      -t "$THREADS" "${extra[@]}" 2>> "$LOG/pileup_${arm}_${tag}.log" \
    | bgzip -@ 4 > "$out"
  if [ ! -s "$out" ]; then echo "[$(stamp)] ⚠️ $arm $tag EMPTY — stopping"; return 1; fi
  tabix -p bed "$out"
  echo "[$(stamp)] --- $arm $tag done, $((( SECONDS-t0 )/60)) min, $(du -h "$out" | cut -f1), $(zcat "$out" | wc -l) rows ---"
}

for arm in NT PA; do
  pileup "$arm" 0.7  gw_thr070 || exit 1
  pileup "$arm" 0.95 gw_thr095 || exit 1
  pileup "$arm" 0.9  feat_thr090 "$FEAT/union_features.bed" || exit 1
  pileup "$arm" 0.98 feat_thr098 "$FEAT/union_features.bed" || exit 1
done

echo "[$(stamp)] === PHASE 3 + 4 COMPLETE — stopping before the contrasts, per PLAN §5 ==="
df -h /data | tail -1
