#!/usr/bin/env bash
# Align the esox basecalls to T2T-CHM13v2.0.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 

set -euo pipefail
P=/data/8oxo_project
OUT=/data/8oxo_project/AdMSC_paper/sessions/2026-08-13/results
REF=$P/refs/T2T-CHM13v2.0.fa
MM2=~/bin/minimap2

echo "[$(date +%H:%M:%S)] minimap2 $($MM2 --version)  samtools $(samtools --version | head -1 | awk '{print $2}')"

for arm in NT PA; do
  echo "[$(date +%H:%M:%S)] ===== $arm ====="
  cat "$P/panel_oxo/esox_$arm"/batch0.fastq "$P/panel_oxo/esox_$arm"/batch1.fastq > "$OUT/esox_${arm}_all.fastq"
  n=$(( $(wc -l < "$OUT/esox_${arm}_all.fastq") / 4 ))
  echo "[$(date +%H:%M:%S)] $arm: $n esox reads"
  "$MM2" -ax map-ont -t 8 "$REF" "$OUT/esox_${arm}_all.fastq" 2> "$OUT/esox_${arm}_mm2.log" \
    | samtools sort -@4 -o "$OUT/esox_${arm}_t2t.bam" -
  samtools index "$OUT/esox_${arm}_t2t.bam"
  echo "[$(date +%H:%M:%S)] $arm: mapped $(samtools view -c -F 0x904 "$OUT/esox_${arm}_t2t.bam") primary alignments of $n reads"
  rm -f "$OUT/esox_${arm}_all.fastq"
done
echo "[$(date +%H:%M:%S)] B2_DONE"
