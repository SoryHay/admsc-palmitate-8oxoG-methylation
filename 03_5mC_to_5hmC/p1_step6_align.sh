#!/usr/bin/env bash
# Align the 5mC/5hmC basecalls to T2T-CHM13v2.0.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
samtools fastq -@4 -T MM,ML,MN "$D/${a}_5mC5hmC.bam" \
  | minimap2 -t 16 -ax map-ont -y --secondary=no "$MMI" - \
  | samtools sort -@4 -m 1500M -o "$D/${a}_aln_T2T.bam" -
samtools index "$D/${a}_aln_T2T.bam"
echo "$a aligned: $(samtools view -c -F 0x904 "$D/${a}_aln_T2T.bam") primary mapped reads" >> "$LOG"
