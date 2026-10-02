#!/usr/bin/env bash
# The same at read quality >= 12 and >= 14 (Table S3).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail; source "$(dirname "$0")/p1_common.sh"; SC=$S/scripts
for a in NT PA; do
  samtools view "$D/${a}_5mC5hmC.bam" | awk -F'\t' '{q=""; for(i=12;i<=NF;i++) if($i ~ /^qs:f:/){q=substr($i,6)}; if(q!="") print $1"\t"q}' > "$D/${a}_qs.tsv"
  echo "$a: $(wc -l < "$D/${a}_qs.tsv") reads with qs; median $(cut -f2 "$D/${a}_qs.tsv" | sort -g | awk '{a[NR]=$1} END{print a[int(NR/2)]}')" >> "$LOG"
  for t in 12 14; do
    awk -v t=$t '$2>=t{print $1}' "$D/${a}_qs.tsv" > "$D/${a}_ids_q${t}.txt"
    samtools view -@4 -b -N "$D/${a}_ids_q${t}.txt" -o "$D/${a}_aln_q${t}.bam" "$D/${a}_aln_T2T.bam"
    samtools index "$D/${a}_aln_q${t}.bam"
    "$MODKIT" pileup "$D/${a}_aln_q${t}.bam" "$D/${a}_cpg_thr07_q${t}.bed" --cpg --ref "$REF" \
      --modified-bases 5mC 5hmC --filter-threshold C:0.7 -t 16 --log-filepath "$D/logs/${a}_modkit_q${t}.log"
    echo "$a Q>=$t: $(wc -l < "$D/${a}_ids_q${t}.txt") reads, pileup $(wc -l < "$D/${a}_cpg_thr07_q${t}.bed") rows" >> "$LOG"
  done
done
for t in 12 14; do
  python3 "$SC/p1_step8_delta.py" "$D/NT_cpg_thr07_q${t}.bed" "$D/PA_cpg_thr07_q${t}.bed" \
    "$S/results/CHECK_P1_qfilter_Q${t}_2026-09-03.md" "CHECK: P-1 recomputed on reads with mean Q >= $t"
done
