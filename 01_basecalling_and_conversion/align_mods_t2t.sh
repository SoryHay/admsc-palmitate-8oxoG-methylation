#!/usr/bin/env bash
# modkit update-tags, samtools fastq -T MM,ML, minimap2 -y to T2T-CHM13v2.0, sort and index.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail
PROJ=/data/8oxo_project
DIR=${1:?sample dir containing pass/ (e.g. methyl/JAR_rep)}
MT=${2:-12}
ST=${3:-4}
MMI=$PROJ/refs/T2T-CHM13v2.0.map-ont.mmi
MODKIT=${MODKIT:-$(ls /data/modkit/latest/*/modkit 2>/dev/null | head -1)}
cd "$PROJ"
[ -d "$DIR/pass" ] || { echo "no $DIR/pass" >&2; exit 2; }
[ -n "$MODKIT" ]   || { echo "modkit not found" >&2; exit 2; }

n=$(ls "$DIR"/pass/*.bam 2>/dev/null | wc -l)
echo ">> [align] $DIR : merging $n pass bams" >&2
samtools cat "$DIR"/pass/*.bam > "$DIR/merged.bam"

echo ">> [align] modkit update-tags --mode implicit (legacy Mm/Ml -> MM/ML)" >&2
"$MODKIT" update-tags --mode implicit "$DIR/merged.bam" "$DIR/upd.bam"
rm -f "$DIR/merged.bam"

echo ">> [align] samtools fastq -T MM,ML | minimap2 -ax map-ont -y -t $MT | sort -@ $ST" >&2
samtools fastq -T MM,ML "$DIR/upd.bam" 2>/dev/null \
  | minimap2 -ax map-ont -y -t "$MT" "$MMI" - 2>"$DIR/minimap2.log" \
  | samtools sort -@ "$ST" -o "$DIR/aln_T2T.bam" -
samtools index "$DIR/aln_T2T.bam"

echo ">> [align] DONE. $(samtools view -c -F 0x904 "$DIR/aln_T2T.bam") primary-mapped reads" >&2
echo ">> [align] mod tags survived? $(samtools view "$DIR/aln_T2T.bam" 2>/dev/null | head -100 | grep -o 'MM:Z:[^[:space:]]*' | grep -oE 'C\+m|A\+a' | sort -u | tr '\n' ' ')" >&2
