#!/usr/bin/env bash
# Convert modkit output to PoreMeth2 input.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
P=/data/8oxo_project/AdMSC_DMR
SRC=$P/readlevel_calls
MAP=/data/8oxo_project/refs/chrom_map.tsv
RESORTER=${RESORTER:-$P/scripts/ModkitResorter.sh}
OUT=${OUT:-$SRC/poremeth}
st(){ date '+%H:%M:%S'; }

ARMS=("$@"); [ ${#ARMS[@]} -eq 0 ] && ARMS=(NT PA)
mkdir -p "$OUT"
[ -s "$RESORTER" ]            || { echo "ERROR: $RESORTER missing"; exit 1; }
[ -s "$P/scripts/ParseModkit.pl" ] || { echo "ERROR: ParseModkit.pl missing"; exit 1; }
[ -s "$MAP" ]                 || { echo "ERROR: $MAP missing"; exit 1; }

for arm in "${ARMS[@]}"; do
  in="$SRC/$arm.calls.tsv"
  ad="$OUT/$arm.modkit_adapted.tsv"
  final="$OUT/$arm.modkit_adapted_sorted.entropy.file.tsv"

  [ -s "$in" ] || { echo "ERROR: $in missing"; exit 1; }
  if [ -s "$final" ]; then echo "[$(st)] $arm — entropy file already present, skipping"; continue; fi

  echo "[$(st)] $arm — adapting $(du -h "$in" | cut -f1) / $(wc -l < "$in") rows to the authors' schema"
  printf 'read_id\t.\tref_position\tchrom\t.\tstrand\t.\t.\t.\t.\tprob_5mC\tmod_code\n' > "$ad"
  LC_ALL=C awk -F'\t' -v OFS='\t' '
      NR==FNR { name[$2]=$1; next }
      ($4 in name) {
        p = ($14=="m") ? $13 : 1-$13
        print $1, ".", $3, name[$4], ".", $6, ".", ".", ".", ".", p, "m"
      }' "$MAP" "$in" >> "$ad"
  echo "[$(st)] $arm — adapted file: $(du -h "$ad" | cut -f1), $(( $(wc -l < "$ad") - 1 )) data rows"
  echo "[$(st)] $arm — contigs present: $(tail -n +2 "$ad" | cut -f4 | sort -u | tr '\n' ' ')"

  echo "[$(st)] $arm — running the authors' ModkitResorter.sh, unmodified"
  sh "$RESORTER" "$ad" 2>&1 | sed "s/^/    /"

  if [ -s "$final" ]; then
    echo "[$(st)] $arm — DONE: $(wc -l < "$final") CpG positions -> $final"
    rm -f "$ad"
  else
    echo "[$(st)] $arm — ⚠️ no entropy file produced; leaving $ad in place for inspection"
  fi
done

echo
echo "=== what the authors' chain kept (columns: chrom pos entropy entropy_cov beta beta_cov) ==="
for arm in "${ARMS[@]}"; do
  f="$OUT/$arm.modkit_adapted_sorted.entropy.file.tsv"
  [ -s "$f" ] || continue
  awk -v a="$arm" -F'\t' '{n++; b+=$6; e+=$4; if($6>=3)g3++; if($6>=5)g5++; if($6>=10)g10++}
    END{printf "%-3s %10d positions   mean beta_cov %.2f   mean entropy_cov %.2f   >=3x %d (%.1f%%)   >=5x %d (%.1f%%)   >=10x %d (%.1f%%)\n",
        a,n,b/n,e/n,g3,100*g3/n,g5,100*g5/n,g10,100*g10/n}' "$f"
done
echo "[$(st)] === ADAPTER + AUTHORS' CHAIN COMPLETE ==="
