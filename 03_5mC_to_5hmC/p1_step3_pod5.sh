#!/usr/bin/env bash
# Convert the subset to pod5.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
pod5 convert fast5 "$D/${a}_fast5"/*.fast5 -o "$D/${a}.pod5" -t 8
echo "$a pod5: $(du -sh "$D/${a}.pod5" | cut -f1), $(pod5 inspect summary "$D/${a}.pod5" 2>/dev/null | grep -iE 'read' | head -1 | tr -s ' ')" >> "$LOG"
