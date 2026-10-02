#!/bin/bash
# Download the T2T-CHM13v2.0 reference.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
cd /data/8oxo_project/refs || exit 1
log(){ echo "[$(date +%H:%M:%S)] $*"; }
export PATH="$HOME/bin:$PATH"

log "=== 1. genome ==="
wget -q -c https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.fa.gz -O T2T-CHM13v2.0.fa.gz \
  && log "genome dl ok ($(du -h T2T-CHM13v2.0.fa.gz|cut -f1))" || { log "GENOME DL FAILED"; }
if [ -s T2T-CHM13v2.0.fa.gz ]; then
  gunzip -kf T2T-CHM13v2.0.fa.gz && log "gunzip ok"
  samtools faidx T2T-CHM13v2.0.fa && log "faidx ok"
  log "chrM check:"; grep -iE 'chrM|MT|CM' T2T-CHM13v2.0.fa.fai | head
  log "=== minimap2 map-ont index ==="; minimap2 -x map-ont -d T2T-CHM13v2.0.map-ont.mmi T2T-CHM13v2.0.fa 2>>fetch.log && log "map-ont mmi ok"
  log "=== minimap2 splice index ==="; minimap2 -x splice -d T2T-CHM13v2.0.splice.mmi T2T-CHM13v2.0.fa 2>>fetch.log && log "splice mmi ok"
fi

log "=== 2. annotation (CAT/Liftoff) ==="
if wget -q -c https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/genes/catLiftOffGenesV1.gff3.gz -O T2T-CHM13v2.0.CAT_liftoff.gff3.gz && [ -s T2T-CHM13v2.0.CAT_liftoff.gff3.gz ]; then
  gunzip -kf T2T-CHM13v2.0.CAT_liftoff.gff3.gz && log "CAT liftoff annotation ok ($(du -h T2T-CHM13v2.0.CAT_liftoff.gff3|cut -f1))"
  echo "ANNOT=CAT_liftoff" > annot_source.txt
else
  log "CAT liftoff URL failed -> trying NCBI RefSeq GCF_009914755.1"
  if wget -q -c 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz' -O T2T-CHM13v2.0.RefSeq.gff.gz && [ -s T2T-CHM13v2.0.RefSeq.gff.gz ]; then
    gunzip -kf T2T-CHM13v2.0.RefSeq.gff.gz && log "NCBI RefSeq annotation ok ($(du -h T2T-CHM13v2.0.RefSeq.gff|cut -f1))"
    echo "ANNOT=NCBI_RefSeq_GCF_009914755.1" > annot_source.txt
  else
    log "BOTH ANNOTATION SOURCES FAILED"; echo "ANNOT=NONE" > annot_source.txt
  fi
fi

log "=== 3. RepeatMasker ==="
if wget -q -c https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/repeatMasker/GCA_009914755.4_rmsk.bb -O T2T-CHM13v2.0.rmsk.bigBed && [ -s T2T-CHM13v2.0.rmsk.bigBed ]; then
  log "rmsk bigBed ok ($(du -h T2T-CHM13v2.0.rmsk.bigBed|cut -f1))"
  if command -v bigBedToBed >/dev/null 2>&1; then
    bigBedToBed T2T-CHM13v2.0.rmsk.bigBed T2T-CHM13v2.0.rmsk.bed && log "rmsk->bed ok"
    grep -iE 'LINE/L1|L1HS|L1PA' T2T-CHM13v2.0.rmsk.bed > T2T-CHM13v2.0.L1.bed && log "L1 subset: $(wc -l < T2T-CHM13v2.0.L1.bed) loci"
  else
    log "bigBedToBed NOT installed - keeping .bb, convert later"
  fi
else
  log "rmsk DL FAILED"
fi
log "=== FETCH DONE ==="
