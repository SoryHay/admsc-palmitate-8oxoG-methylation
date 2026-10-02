#!/usr/bin/env bash
# esox 8-oxo-dG calling: resquiggle, bonito and remora on the dorado guide (Docker).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
F5DIR=${1:?fast5 dir}; OUT=${2:?output dir}; shift 2 || true
GPU=0; SUP=/data/dorado/models_r9_v33/dna_r9.4.1_e8_sup@v3.3; IMG=${ESOX_IMG:-esox:latest}
while [ $# -gt 0 ]; do case "$1" in --gpu) GPU=$2; shift 2;; --sup) SUP=$2; shift 2;; *) echo "unknown $1">&2; exit 2;; esac; done
F5DIR=$(realpath "$F5DIR"); OUT=$(realpath -m "$OUT"); mkdir -p "$OUT"/guide "$OUT"/esox "$OUT"/modcall

for f5 in "$F5DIR"/*.fast5; do
  base=$(basename "$f5" .fast5)
  echo ">> [esox_run] SUP guide: $base" >&2
  "$HERE/dorado_docker.sh" "$SUP" "$f5" "$OUT/guide/${base}.fastq" --gpu "$GPU" ${DORADO_EXTRA:-}
done

edock(){ docker run --rm --name "$1" --gpus "device=$GPU" -v /data/8oxo_project:/data/8oxo_project \
    -v "$F5DIR":"$F5DIR":ro -v "$OUT":"$OUT" -w /opt/esox "$IMG" \
    conda run -n esox_env python3 "${@:2}"; }

echo ">> [esox_run] esox basecall (bonito, docker)" >&2
edock esoxbc_$$ scripts/basecall.py --fast5-path "$F5DIR" --fastq-path "$OUT/guide" \
  --output-path "$OUT/esox" --model-file static/models/bonito.pt --device cuda:0

echo ">> [esox_run] esox modcall (remora, docker)" >&2
edock esoxmc_$$ scripts/modcall.py --input-path "$OUT/esox" --output-path "$OUT/modcall" \
  --model-file static/models/remora.pt --device cuda:0

echo ">> [esox_run] DONE. modcall rows: $(cat "$OUT"/modcall/*.txt 2>/dev/null | wc -l)  (image $(docker inspect --format='{{.Id}}' "$IMG" 2>/dev/null|cut -c1-19))" >&2
