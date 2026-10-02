#!/usr/bin/env bash
# CpG-only bedMethyl from the context table.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
B=/data/8oxo_project/AdMSC_DMR/trackB_rerio_allctx
for arm in NT PA; do
  o=$B/$arm/cpg_only_thr070.bed.gz
  if [ -s "$o.tbi" ]; then echo "$arm already present"; continue; fi
  zcat "$B/$arm/ctx_thr070.cpg.tsv.gz" | tail -n +2 | awk 'BEGIN{OFS="\t"}
    {nv=$5; nm=$6; nc=nv-nm; f=(nv>0)?100*nm/nv:0;
     print $1,$2,$3,"m",nv,$4,$2,$3,"255,0,0",nv,sprintf("%.2f",f),nm,nc,0,0,0,0,0}' \
    | bgzip -@ 4 > "$o"
  tabix -p bed "$o"
  echo "$arm: $(zcat "$o" | wc -l) CpG positions"
done
