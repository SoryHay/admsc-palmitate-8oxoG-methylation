#!/usr/bin/env bash
# Extract the subset reads from fast5_pass.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
fast5_subset -i "${FAST5[$a]}" -s "$D/${a}_fast5" -l "$D/${a}_ids_10pct.txt" -t 8 -n 4000
n=$(grep -oE '[0-9]+ reads extracted' "$D/logs/step2_fast5_${a}.log" 2>/dev/null | tail -1 || true)
echo "$a subset: ${n:-see log} / $(wc -l < "$D/${a}_ids_10pct.txt") requested, $(du -sh "$D/${a}_fast5" | cut -f1)" >> "$LOG"
