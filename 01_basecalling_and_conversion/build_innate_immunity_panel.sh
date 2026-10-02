#!/usr/bin/env bash
# Innate-immunity gene panel at exon and span level.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail
PROJ=/data/8oxo_project
FEAT=$PROJ/AdMSC_DMR/features
GFF=$FEAT/.t2t_refseq.cp.gff
LIST=$FEAT/innate_immunity_genes.tsv

[ -s "$GFF" ] || { echo "run build_features.sh first (translated GFF missing)"; exit 1; }

if [ ! -s "$LIST" ]; then
cat > "$LIST" <<'TSV'
TLR1	TLR
TLR2	TLR
TLR3	TLR
TLR4	TLR
TLR5	TLR
TLR6	TLR
TLR7	TLR
TLR8	TLR
TLR9	TLR
TLR10	TLR
CD14	TLR
LY96	TLR
MYD88	TLR_adaptor
TIRAP	TLR_adaptor
TICAM1	TLR_adaptor
TICAM2	TLR_adaptor
IRAK1	TLR_adaptor
IRAK2	TLR_adaptor
IRAK3	TLR_adaptor
IRAK4	TLR_adaptor
TRAF3	TLR_adaptor
TRAF6	TLR_adaptor
NLRP1	NLR
NLRP2	NLR
NLRP3	NLR
NLRP4	NLR
NLRP5	NLR
NLRP6	NLR
NLRP7	NLR
NLRP8	NLR
NLRP9	NLR
NLRP10	NLR
NLRP11	NLR
NLRP12	NLR
NLRP13	NLR
NLRP14	NLR
NLRC3	NLR
NLRC4	NLR
NLRC5	NLR
NLRX1	NLR
NOD1	NLR
NOD2	NLR
NAIP	NLR
CIITA	NLR
DDX58	RLR
IFIH1	RLR
DHX58	RLR
MAVS	RLR
CGAS	cGAS_STING
STING1	cGAS_STING
TBK1	cGAS_STING
IKBKE	cGAS_STING
IRF3	cGAS_STING
IRF7	cGAS_STING
AIM2	inflammasome
IFI16	inflammasome
PYCARD	inflammasome
CASP1	inflammasome
CASP4	inflammasome
CASP5	inflammasome
CASP8	inflammasome
GSDMD	inflammasome
GSDME	inflammasome
NEK7	inflammasome
MEFV	inflammasome
PSTPIP1	inflammasome
IL1A	cytokine_IL1
IL1B	cytokine_IL1
IL1RN	cytokine_IL1
IL18	cytokine_IL1
IL18BP	cytokine_IL1
IL33	cytokine_IL1
IL36A	cytokine_IL1
IL36B	cytokine_IL1
IL36G	cytokine_IL1
IL37	cytokine_IL1
IL1F10	cytokine_IL1
IL2	cytokine
IL4	cytokine
IL6	cytokine
IL7	cytokine
IL10	cytokine
IL11	cytokine
IL12A	cytokine
IL12B	cytokine
IL13	cytokine
IL15	cytokine
IL17A	cytokine
IL17F	cytokine
IL21	cytokine
IL22	cytokine
IL23A	cytokine
IL27	cytokine
IL34	cytokine
TNF	TNF_superfamily
LTA	TNF_superfamily
LTB	TNF_superfamily
TNFSF10	TNF_superfamily
TNFSF11	TNF_superfamily
CD40LG	TNF_superfamily
FASLG	TNF_superfamily
TNFRSF1A	TNF_superfamily
TNFRSF1B	TNF_superfamily
TNFAIP3	NFkB
NFKB1	NFkB
NFKB2	NFkB
RELA	NFkB
RELB	NFkB
REL	NFkB
NFKBIA	NFkB
CHUK	NFkB
IKBKB	NFkB
IKBKG	NFkB
IFNA1	interferon
IFNA2	interferon
IFNB1	interferon
IFNG	interferon
IFNL1	interferon
IFNL2	interferon
IFNL3	interferon
IFNAR1	interferon
IFNAR2	interferon
IFNGR1	interferon
IFNGR2	interferon
IFNLR1	interferon
ISG15	ISG
IFIT1	ISG
IFIT2	ISG
IFIT3	ISG
IFITM1	ISG
IFITM3	ISG
MX1	ISG
MX2	ISG
OAS1	ISG
OAS2	ISG
OAS3	ISG
OASL	ISG
RSAD2	ISG
BST2	ISG
ZBP1	ISG
STAT1	ISG
STAT2	ISG
IRF1	ISG
IRF9	ISG
CCL2	chemokine
CCL3	chemokine
CCL4	chemokine
CCL5	chemokine
CCL20	chemokine
CXCL1	chemokine
CXCL2	chemokine
CXCL8	chemokine
CXCL9	chemokine
CXCL10	chemokine
CXCL11	chemokine
CXCL12	chemokine
CCR2	chemokine
CCR5	chemokine
CXCR4	chemokine
TGFB1	growth_factor
TGFB2	growth_factor
TGFB3	growth_factor
CSF1	growth_factor
CSF2	growth_factor
CSF3	growth_factor
CD36	lipid_sensing
MSR1	lipid_sensing
CLEC7A	lipid_sensing
CLEC4E	lipid_sensing
MRC1	lipid_sensing
TREM1	lipid_sensing
TREM2	lipid_sensing
FFAR1	lipid_sensing
FFAR4	lipid_sensing
TSV
fi

NWANT=$(grep -vc '^#' "$LIST")
echo "curated symbols: $NWANT"

awk -F'\t' 'NR==FNR{ if($0 !~ /^#/){split($0,a,"\t"); fam[a[1]]=a[2]} next }
     ($4 in fam){ print $1,$2,$3,$4,fam[$4],$6 }' OFS='\t' \
     "$LIST" "$FEAT/genebodies.bed" | sort -k1,1 -k2,2n > "$FEAT/innate_genebodies.bed"

awk 'BEGIN{OFS="\t"}{ if($6=="-") tss=$3; else tss=$2; s=tss-1000; if(s<0)s=0;
     print $1,s,tss+1000,$4"_prom",$5,$6 }' "$FEAT/innate_genebodies.bed" \
  | sort -k1,1 -k2,2n > "$FEAT/innate_promoters.bed"

awk -F'\t' 'NR==FNR{ if($0 !~ /^#/){split($0,a,"\t"); fam[a[1]]=a[2]} next }
     { n=split($4,g,","); for(i=1;i<=n;i++) if (g[i] in fam){ print $1,$2,$3,g[i],fam[g[i]]; break } }' \
     OFS='\t' "$LIST" "$FEAT/exons_merged.bed" | sort -k1,1 -k2,2n > "$FEAT/innate_exons.bed"

FOUND=$(cut -f4 "$FEAT/innate_genebodies.bed" | sort -u | wc -l)
echo "found in annotation: $FOUND / $NWANT"
echo
echo "MISSING (reported, not dropped silently):"
comm -23 <(grep -v '^#' "$LIST" | cut -f1 | sort -u) \
         <(cut -f4 "$FEAT/innate_genebodies.bed" | sort -u) | tr '\n' ' '
echo; echo
echo "per family (genes / exon intervals):"
awk -F'\t' '{g[$5"\t"$4]=1; e[$5]++} END{for(k in g){split(k,a,"\t"); n[a[1]]++}
     for(f in e) printf "  %-18s %3d genes  %5d exons\n", f, n[f], e[f]}' \
     "$FEAT/innate_exons.bed" | sort
echo
printf '  %-26s %6d intervals\n' innate_genebodies "$(wc -l < "$FEAT/innate_genebodies.bed")" \
    innate_promoters "$(wc -l < "$FEAT/innate_promoters.bed")" \
    innate_exons "$(wc -l < "$FEAT/innate_exons.bed")"
