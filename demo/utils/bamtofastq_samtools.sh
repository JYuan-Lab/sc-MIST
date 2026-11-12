#!/usr/bin/env bash
# Minimal drop-in replacement for 10x bamtofastq, used only by the demo so
# that the demo does not require a Cell Ranger installation.
set -euo pipefail

THREADS=1
BAM=""
OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --nthreads) THREADS="$2"; shift 2 ;;
    --nthreads=*) THREADS="${1#*=}"; shift ;;
    -*) shift ;;
    *) if [[ -z "$BAM" ]]; then BAM="$1"; elif [[ -z "$OUTDIR" ]]; then OUTDIR="$1"; fi; shift ;;
  esac
done

[[ -s "$BAM" && -n "$OUTDIR" ]] || { echo "usage: bamtofastq_samtools.sh --nthreads N in.bam outdir" >&2; exit 1; }
mkdir -p "$OUTDIR"

samtools sort -n -@ "$THREADS" -o "$OUTDIR/name_sorted.bam" "$BAM"
samtools fastq -@ "$THREADS" -n \
  -1 "$OUTDIR/demo_R1.fastq" \
  -2 "$OUTDIR/demo_R2.fastq" \
  -0 /dev/null -s /dev/null \
  "$OUTDIR/name_sorted.bam"
gzip -f "$OUTDIR/demo_R1.fastq" "$OUTDIR/demo_R2.fastq"
