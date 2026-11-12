#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Refresh read-to-taxid assignments on host-filtered reads. This is an
internal helper for step 05: downstream barcode/UMI-aware steps need a fresh
read-to-taxid table after HLA/rRNA filtering.

Required:
  --sample        Sample ID
  --r1            Clean R1 FASTQ
  --r2            Clean R2 FASTQ
  --outdir        Output directory (step 05 sample directory)
  --db            Kraken2 database directory

Optional:
  --kraken2       kraken2 or Kraken2Uniq (default: kraken2)
  --threads       Threads (default: 8)
  --rscript       Rscript executable (default: Rscript)
  --script-dir    Directory containing the 02/03 helper scripts (default: this script's dir)
EOF
}

SAMPLE=""
R1=""
R2=""
OUTDIR=""
DB=""
KRAKEN2="kraken2"
THREADS=8
RSCRIPT="Rscript"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --db) DB="$2"; shift 2 ;;
    --kraken2) KRAKEN2="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLE" && -n "$R1" && -n "$R2" && -n "$OUTDIR" && -n "$DB" ]] || {
  usage >&2
  exit 1
}

mkdir -p "$OUTDIR" "$OUTDIR/microbial"

bash "$SCRIPT_DIR/02_kraken_classify.sh" \
  --sample "$SAMPLE" --r1 "$R1" --r2 "$R2" \
  --outdir "$OUTDIR" --db "$DB" --kraken2 "$KRAKEN2" --threads "$THREADS"

"$RSCRIPT" "$SCRIPT_DIR/03_extract_microbial_reads.R" \
  --sample "$SAMPLE" \
  --kraken-report "$OUTDIR/${SAMPLE}.kraken.report.txt" \
  --mpa-report "$OUTDIR/${SAMPLE}.kraken.report.mpa.txt" \
  --kraken-output "$OUTDIR/${SAMPLE}.kraken.output.txt" \
  --fq1 "$OUTDIR/${SAMPLE}_1.fq" \
  --fq2 "$OUTDIR/${SAMPLE}_2.fq" \
  --outdir "$OUTDIR/microbial"

echo "[05] done: $SAMPLE"
