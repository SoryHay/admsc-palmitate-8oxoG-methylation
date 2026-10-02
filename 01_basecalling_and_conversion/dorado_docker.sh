#!/usr/bin/env bash
# dorado 0.9.6 SUP guide basecall for esox (Docker); FASTQ headers reduced to the read ID.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -euo pipefail

MODEL=${1:?model path (under /data/dorado)}; IN=${2:?input fast5/pod5/dir}; OUT=${3:?output file}; shift 3 || true
DORADO_BIN=/data/dorado/dorado-0.9.6-linux-x64/bin/dorado
GPU=0
IMG=${DORADO_IMG:-esox:latest}
EMIT=fastq
EXTRA=()
while [ $# -gt 0 ]; do case "$1" in
  --bam) EMIT=bam; shift;;
  --dorado) DORADO_BIN=$2; shift 2;;
  --gpu) GPU=$2; shift 2;;
  *) EXTRA+=("$1"); shift;;
esac; done

IN=$(realpath "$IN"); OUTDIR=$(dirname "$(realpath -m "$OUT")"); mkdir -p "$OUTDIR"
NAME="dorado_$(basename "$OUT" | tr -c 'A-Za-z0-9_.-' '_')_$$"
echo "$NAME" > "${OUT}.container"          # so a hung run is killable: docker kill $(cat OUT.container)
echo ">> dorado in docker  [$NAME]  GPU device=$GPU  model=$(basename "$MODEL")" >&2
echo "   kill if hung:  docker kill $NAME" >&2

run_docker(){ docker run --rm --name "$NAME" --gpus "device=$GPU" \
    -v /data/dorado:/data/dorado:ro -v "$IN":"$IN":ro -v "$OUTDIR":"$OUTDIR" \
    --entrypoint "$DORADO_BIN" "$IMG" \
    basecaller "$MODEL" "$IN" --device cuda:0 "$@"; }

if [ "$EMIT" = bam ]; then
  run_docker "${EXTRA[@]}" > "$OUT" 2> "${OUT}.log"
else
  run_docker --emit-fastq "${EXTRA[@]}" 2> "${OUT}.log" \
    | awk 'NR%4==1{print $1; next} {print}' > "$OUT"
fi
rm -f "${OUT}.container"
echo ">> done -> $OUT ($(grep -oE 'Simplex reads basecalled: [0-9]+' "${OUT}.log" | tail -1))" >&2
