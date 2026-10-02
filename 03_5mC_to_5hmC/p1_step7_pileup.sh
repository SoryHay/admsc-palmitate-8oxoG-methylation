#!/usr/bin/env bash
# modkit pileup --cpg with 5mC and 5hmC at threshold C:0.7.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; a=${1:?NT|PA}
"$MODKIT" pileup "$D/${a}_aln_T2T.bam" "$D/${a}_cpg_thr07.bed" --cpg --ref "$REF" \
  --modified-bases 5mC 5hmC --filter-threshold C:0.7 -t 16 --log-filepath "$D/logs/${a}_modkit.log"
echo "$a pileup: $(wc -l < "$D/${a}_cpg_thr07.bed") rows" >> "$LOG"
