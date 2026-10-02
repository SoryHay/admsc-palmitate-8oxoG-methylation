# Shared settings of the 5mC/5hmC pipeline (10 % read subset).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
P=/data/8oxo_project
S=$P/AdMSC_paper/sessions/2026-09-03
D=$S/p1_5hmC
LOG=$D/RUN_LOG.md
MODEL=/data/dorado/models_r9_v33/dna_r9.4.1_e8_sup@v3.3
MOD=/data/dorado/models_r9_v33/dna_r9.4.1_e8_sup@v3.3_5mCG_5hmCG@v0
MODKIT=/data/modkit/latest/dist_modkit_v0.6.4_cd85862/modkit
REF=$P/refs/T2T-CHM13v2.0.fa
MMI=$P/refs/T2T-CHM13v2.0.map-ont.mmi
NAS=<DATA_DISK>/ONT_Seq_for_Analysis/AdMSC_nt_PA_1
declare -A FAST5=( [NT]=$NAS/NT-AdMSCs/fast5_pass [PA]=$NAS/250uM24h-AdMSCs/fast5_pass )
