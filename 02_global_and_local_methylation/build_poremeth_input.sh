#!/usr/bin/env bash
# Build the PoreMeth2 input for both arms.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
P=/data/8oxo_project/AdMSC_DMR
T=$P/readlevel_calls
R=$P/results
MINPROB=0.70
st(){ date '+%H:%M:%S'; }

mkdir -p "$R"

for arm in NT PA; do
  out="$T/$arm.cpg_beta.tsv"
  if [ -s "$out" ]; then echo "[$(st)] $arm.cpg_beta.tsv exists, skipping"; continue; fi
  echo "[$(st)] $arm — filter (CpG, call_prob >= $MINPROB), sort, aggregate per position"
  LC_ALL=C awk -F'\t' -v mp="$MINPROB" '
      ($13+0) >= mp && substr($16,3,1)=="c" && substr($16,4,1)=="g" {
          print $4"\t"$3"\t"($14=="m"?1:0)
      }' "$T/$arm.calls.tsv" \
    | LC_ALL=C sort -k1,1 -k2,2n -S 4G --parallel=8 -T "$T" \
    | LC_ALL=C awk -F'\t' '
        BEGIN{OFS="\t"; c=""; p=-1; m=0; n=0}
        {
          if ($1!=c || $2!=p) {
            if (n>0) print c,p,m,n
            c=$1; p=$2; m=0; n=0
          }
          m+=$3; n++
        }
        END{ if (n>0) print c,p,m,n }' > "$out"
  echo "[$(st)] $arm — $(wc -l < "$out") CpG positions retained"
done

echo "[$(st)] joining arms and measuring FW feasibility"
LC_ALL=C join -t $'\t' -j 1 \
  <(LC_ALL=C awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$T/NT.cpg_beta.tsv") \
  <(LC_ALL=C awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$T/PA.cpg_beta.tsv") \
  > "$T/both.cpg_beta.tsv"
echo "[$(st)] CpG positions covered in BOTH arms: $(wc -l < "$T/both.cpg_beta.tsv")"

python3 - "$T/both.cpg_beta.tsv" "$R/poremeth_input_feasibility.txt" <<'PY'
import sys
from collections import Counter

src, out = sys.argv[1], sys.argv[2]
COVS = [1, 2, 3, 4, 5, 10]
GAPS = [500, 1000, 5000]

rows = []
with open(src) as f:
    for line in f:
        k, ntm, ntc, pam, pac = line.rstrip('\n').split('\t')
        c, p = k.rsplit(':', 1)
        rows.append((c, int(p), int(ntc), int(pac)))
rows.sort(key=lambda r: (r[0], r[1]))

lines = []
w = lines.append
w("PoreMeth2 INPUT FEASIBILITY — can FW >= 3 be satisfied?")
w(f"source: {src}")
w(f"operating point: call_prob >= 0.70 (locked, PLAN.md §21)")
w(f"CpG positions covered in both arms: {len(rows):,}")
w("")
w("A CpG is KEPT if coverage >= t in BOTH arms. A RUN is a maximal stretch of kept CpGs whose")
w("consecutive genomic gaps are all <= gap. FW = 3 needs runs of length >= 3.")
w("")
hdr = f"{'cov>=t':>7}{'kept CpG':>14}{'% of both':>11}"
for g in GAPS:
    hdr += f"{'runs>=3 (g'+str(g)+')':>20}{'CpG in them':>14}"
w(hdr)

for t in COVS:
    kept = [r for r in rows if r[2] >= t and r[3] >= t]
    line = f"{t:>7}{len(kept):>14,}{100*len(kept)/max(1,len(rows)):>10.2f}%"
    for g in GAPS:
        runs, cur = [], []
        prevc, prevp = None, None
        for c, p, _, _ in kept:
            if prevc == c and p - prevp <= g:
                cur.append(p)
            else:
                if len(cur) >= 3: runs.append(len(cur))
                cur = [p]
            prevc, prevp = c, p
        if len(cur) >= 3: runs.append(len(cur))
        line += f"{len(runs):>20,}{sum(runs):>14,}"
    w(line)

w("")
w("Read the table as: at threshold t, DMR calling is POSSIBLE only inside the runs counted here.")
w("A high 'kept CpG' with few runs means the survivors are scattered and the segmentation would")
w("be joining positions across unmeasured sequence.")

open(out, 'w').write("\n".join(lines) + "\n")
print("\n".join(lines))
PY
echo "[$(st)] === DONE -> $R/poremeth_input_feasibility.txt ==="
