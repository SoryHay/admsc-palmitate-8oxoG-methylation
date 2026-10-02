#!/usr/bin/env bash
# The same by decile of read quality (Table S4).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; SC=$S/scripts
for a in NT PA; do
  "$MODKIT" extract calls "$D/${a}_aln_T2T.bam" "$D/${a}_calls_cpg.tsv" --cpg --reference "$REF" \
    --filter-threshold C:0.7 --pass-only --mapped-only -t 16 --log-filepath "$D/logs/${a}_extract.log"
  echo "$a extract calls: $(wc -l < "$D/${a}_calls_cpg.tsv") rows" >> "$LOG"
done
python3 "$SC/p1_check_qdecile.py" "$D" "$S/results/CHECK_P1_qdecile_2026-09-03.md"
