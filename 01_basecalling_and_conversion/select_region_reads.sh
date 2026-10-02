#!/usr/bin/env bash
# Read IDs whose primary alignment overlaps a BED region set.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail
BED=${1:?region bed (matching the bam reference)}; OUT=${2:?output read-ids file}; shift 2 || true
[ $# -ge 1 ] || { echo "give a bam dir or bam files" >&2; exit 2; }

BAMS=()
for x in "$@"; do
  if [ -d "$x" ]; then while IFS= read -r b; do BAMS+=("$b"); done < <(find "$x" -maxdepth 1 -iname '*.bam');
  else BAMS+=("$x"); fi
done
echo ">> [select_region_reads] ${#BAMS[@]} bams vs $(wc -l < "$BED") regions -> $OUT" >&2

: > "$OUT.tmp"
for b in "${BAMS[@]}"; do
  samtools view -L "$BED" "$b" 2>/dev/null | cut -f1
done | sort -u > "$OUT"
rm -f "$OUT.tmp"
echo ">> [select_region_reads] $(wc -l < "$OUT") unique reads overlap the regions" >&2
