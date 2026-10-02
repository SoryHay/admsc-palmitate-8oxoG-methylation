#!/usr/bin/env bash
# Choose the dorado batch size so that peak GPU memory stays below 80 %.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
"$P/scripts/dorado_docker.sh" "$MODEL" "$D/${a}.pod5" "$D/bench_${a}.bam" --bam --gpu 0 \
  --modified-bases-models "$MOD" --verbose -b 0 -n 200
LOGF=$D/bench_${a}.bam.log
LIMIT=$(grep -m1 -oE 'memory limit [0-9.]+GB' "$LOGF" | grep -oE '[0-9.]+' || true)
BATCH=$(grep -m1 -oE 'using chunk size [0-9]+, batch size [0-9]+' "$LOGF" | grep -oE '[0-9]+$' || true)
MODMEM=$(grep -m1 -oE 'Model memory [0-9.]+GB' "$LOGF" | grep -oE '[0-9.]+' || true)
DECMEM=$(grep -m1 -oE 'Decode memory [0-9.]+GB' "$LOGF" | grep -oE '[0-9.]+' || true)
if [ -z "$LIMIT" ] || [ -z "$BATCH" ] || [ -z "$MODMEM" ] || [ -z "$DECMEM" ]; then
  echo "bench parse FAILED: limit='$LIMIT' batch='$BATCH' model='$MODMEM' decode='$DECMEM' — see $LOGF" | tee -a "$LOG"; exit 1; fi
awk -v lim="$LIMIT" -v b="$BATCH" -v mm="$MODMEM" -v dm="$DECMEM" 'BEGIN{used=mm+dm; target=0.80*lim;
  if (used>target) { rb=int(b*target/used/64)*64; if(rb<64)rb=64 } else { rb=b }; print rb}' > "$D/batch_pin.txt"
echo "bench ($a, 200 reads, docker GPU0): card limit ${LIMIT} GB, auto batch ${BATCH}, model ${MODMEM} GB + decode ${DECMEM} GB -> pin -b $(cat "$D/batch_pin.txt")" >> "$LOG"
rm -f "$D/bench_${a}.bam"
