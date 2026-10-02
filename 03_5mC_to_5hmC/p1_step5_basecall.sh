#!/usr/bin/env bash
# dorado 0.9.6 with dna_r9.4.1_e8_sup@v3.3 and 5mCG_5hmCG@v0 (Docker, GPU).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
B=$(cat "$D/batch_pin.txt"); [ -n "$B" ] || { echo "no batch pin" | tee -a "$LOG"; exit 1; }
"$P/scripts/dorado_docker.sh" "$MODEL" "$D/${a}.pod5" "$D/${a}_5mC5hmC.bam" --bam --gpu 0 \
  --modified-bases-models "$MOD" -b "$B"
echo "$a basecall: $(grep -oE 'Simplex reads basecalled: [0-9]+' "$D/${a}_5mC5hmC.bam.log" | tail -1), $(grep -oE 'Finished in \(ms\): [0-9]+' "$D/${a}_5mC5hmC.bam.log" | tail -1), -b $B" >> "$LOG"
