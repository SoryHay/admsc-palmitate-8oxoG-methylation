#!/usr/bin/env bash
# Region annotation on T2T: CpG islands, shores, shelves, promoters, exons, gene bodies.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail

PROJ=/data/8oxo_project
FEAT=$PROJ/AdMSC_DMR/features
GFF=$PROJ/refs/T2T_RefSeq.gff
MAP=$PROJ/refs/chrom_map.tsv
FAI=$PROJ/refs/T2T-CHM13v2.0.fa.fai
ISL=$PROJ/refs/T2T-CHM13v2.0.cpgIslands.strict.bed
mkdir -p "$FEAT"

cut -f1,2 "$FAI" | sort -k1,1 > "$FEAT/t2t.genome"

echo "== islands (existing, not recomputed) =="
ln -sf "$ISL" "$FEAT/islands.bed"
echo "   $(wc -l < "$ISL") islands"

echo "== shores: 2 kb flanks, island body removed =="
sort -k1,1 -k2,2n "$ISL" | cut -f1-4 > "$FEAT/.isl.sorted.bed"
bedtools flank -i "$FEAT/.isl.sorted.bed" -g "$FEAT/t2t.genome" -b 2000 \
  | bedtools subtract -a - -b "$FEAT/.isl.sorted.bed" \
  | awk 'BEGIN{OFS="\t"} $3-$2>=200 {print $1,$2,$3,$4"_shore"}' \
  | sort -k1,1 -k2,2n > "$FEAT/shores.bed"
echo "   $(wc -l < "$FEAT/shores.bed") shore segments"

echo "== translate RefSeq NC_* -> CP* and extract gene/exon features =="
awk -F'\t' 'NR==FNR{m[$3]=$2; next} /^#/{next} ($1 in m){$1=m[$1]; print}' OFS='\t' \
    "$MAP" "$GFF" > "$FEAT/.t2t_refseq.cp.gff"
echo "   GFF rows kept: $(wc -l < "$FEAT/.t2t_refseq.cp.gff")  (of $(grep -vc '^#' "$GFF"))"

awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $9~/gene_biotype=protein_coding/ {
      name="."; if (match($9,/Name=[^;]+/)) name=substr($9,RSTART+5,RLENGTH-5);
      print $1,$4-1,$5,name,".",$7 }' "$FEAT/.t2t_refseq.cp.gff" \
  | sort -k1,1 -k2,2n > "$FEAT/genebodies.bed"
echo "   $(wc -l < "$FEAT/genebodies.bed") protein-coding gene bodies"

awk 'BEGIN{OFS="\t"} { if ($6=="-") {tss=$3} else {tss=$2};
      s=tss-1000; if (s<0) s=0; print $1,s,tss+1000,$4"_prom",".",$6 }' \
    "$FEAT/genebodies.bed" | sort -k1,1 -k2,2n > "$FEAT/promoters_tss1kb.bed"
echo "   $(wc -l < "$FEAT/promoters_tss1kb.bed") promoters (TSS +-1 kb)"

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {
      name="."; if (match($9,/gene=[^;]+/)) name=substr($9,RSTART+5,RLENGTH-5);
      print $1,$4-1,$5,name }' "$FEAT/.t2t_refseq.cp.gff" \
  | sort -k1,1 -k2,2n -k4,4 \
  | bedtools merge -i - -c 4 -o distinct > "$FEAT/exons_merged.bed"
echo "   $(wc -l < "$FEAT/exons_merged.bed") merged exon intervals"

echo "== 82-gene panel, at EXON level (never span) =="
if [ -s "$PROJ/refs/panel_t2t.bed" ]; then
   grep -v "^#" "$PROJ/refs/panel_t2t.bed" | cut -f5 | sort -u > "$FEAT/.panel_genes.txt"
  awk -F'\t' 'NR==FNR{g[$1]; next} { n=split($4,a,","); for(i=1;i<=n;i++) if (a[i] in g) {print; break} }' \
      OFS='\t' "$FEAT/.panel_genes.txt" "$FEAT/exons_merged.bed" > "$FEAT/panel_exons.bed"
  echo "   $(wc -l < "$FEAT/panel_exons.bed") panel exon intervals from $(wc -l < "$FEAT/.panel_genes.txt") gene names"
fi

rm -f "$FEAT/.isl.sorted.bed"
echo
echo "== summary =="
for f in islands shores promoters_tss1kb genebodies exons_merged panel_exons; do
  [ -e "$FEAT/$f.bed" ] || continue
  printf '  %-20s %8d intervals  %12d bp\n' "$f" \
     "$(wc -l < "$FEAT/$f.bed")" \
     "$(awk '{s+=$3-$2} END{print s+0}' "$FEAT/$f.bed")"
done
