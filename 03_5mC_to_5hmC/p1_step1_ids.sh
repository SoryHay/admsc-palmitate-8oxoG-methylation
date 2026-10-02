#!/usr/bin/env bash
# Seeded 10 % random subset of pass read IDs per arm.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail
P=/data/8oxo_project; D=$P/AdMSC_paper/sessions/2026-09-03/p1_5hmC
for a in NT PA; do
  samtools view -@4 -F 0x900 $P/AdMSC_DMR/trackB_rerio_allctx/$a/aln_T2T.bam | cut -f1 | LC_ALL=C sort -u -S 2G > $D/${a}_all_ids.txt
  python3 - "$D/${a}_all_ids.txt" "$D/${a}_ids_10pct.txt" 20260903 <<'PY'
import random, sys
ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
k = len(ids) // 10
sub = random.Random(int(sys.argv[3])).sample(ids, k)
open(sys.argv[2], 'w').write('\n'.join(sub) + '\n')
print(f"{sys.argv[1]}: {len(ids)} ids -> {k} sampled (seed {sys.argv[3]})")
PY
done
