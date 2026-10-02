#!/usr/bin/env bash
# Guppy 5.0.16 + Rerio all-context basecalling of both arms (Docker, GPU).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail
PROJ=/data/8oxo_project
ONT=<DATA_DISK>/ONT_Seq_for_Analysis/AdMSC_nt_PA_1
RERIO=~/rerio
CFG=res_dna_r941_min_modbases-all-context_v001.cfg
IMG=genomicpariscentre/guppy-gpu:5.0.16
GPU=${1:-0}
OUT=$PROJ/AdMSC_DMR/trackB_rerio_allctx
LOG=$PROJ/AdMSC_DMR/logs
mkdir -p "$OUT" "$LOG"

stamp(){ date '+%Y-%m-%d %H:%M:%S'; }

run_arm(){                       # $1 = NT|PA   $2 = fast5 dir
  local arm=$1 f5=$2 o="$OUT/$1"
  mkdir -p "$o"
  echo "[$(stamp)] === $arm START — $(ls "$f5"/*.fast5 | wc -l) fast5, $(du -sh "$f5" | cut -f1) ==="
  echo "[$(stamp)] kill with: docker kill guppy_AdMSC_$arm"
  local t0=$SECONDS
  docker run --rm --name "guppy_AdMSC_$arm" --gpus "device=$GPU" \
    -v "$RERIO":/rerio:ro -v "$f5":/in:ro -v "$o":/out "$IMG" \
    guppy_basecaller -i /in -s /out -c "$CFG" -d /rerio/basecall_models \
    -x cuda:0 --bam_out > "$LOG/basecall_$arm.guppy.log" 2>&1
  local rc=$? dt=$((SECONDS-t0))
  local reads pass fail tags
  reads=$(awk 'NR>1{n++} END{print n+0}' "$o/sequencing_summary.txt" 2>/dev/null)
  pass=$(ls "$o/pass"/*.bam 2>/dev/null | wc -l)
  fail=$(ls "$o/fail"/*.bam 2>/dev/null | wc -l)
  tags=$(samtools view "$(ls "$o/pass"/*.bam 2>/dev/null | head -1)" 2>/dev/null | head -2000 | grep -c 'Mm:Z:')
  echo "[$(stamp)] === $arm DONE — rc=$rc  wall=$((dt/60)) min  reads=$reads  pass_bams=$pass fail_bams=$fail  tags=$tags/2000 ==="
  echo "$arm rc=$rc wall_min=$((dt/60)) reads=$reads pass_bams=$pass tags=$tags" >> "$LOG/basecall_summary.txt"
  [ "$rc" -eq 0 ] && [ "${reads:-0}" -gt 0 ] && [ "$tags" -gt 1900 ]
}

echo "[$(stamp)] Phase 2 launched on GPU device=$GPU"
if run_arm NT "$ONT/NT-AdMSCs/fast5_pass"; then
  echo "[$(stamp)] NT PASSED the acceptance criteria -> starting PA"
  run_arm PA "$ONT/250uM24h-AdMSCs/fast5_pass" \
    && echo "[$(stamp)] BOTH ARMS COMPLETE" \
    || echo "[$(stamp)] ⚠️ PA FAILED acceptance — stop and review"
else
  echo "[$(stamp)] ⚠️ NT FAILED acceptance — PA NOT started, per plan §2"
fi
