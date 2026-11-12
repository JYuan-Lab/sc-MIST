#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Step 05: remove HLA/MHC and human rRNA reads, then refresh the
read-to-taxid assignments on the reads that survived.

Usage:
  bash scripts/05_filter_host_and_refresh.sh \
    --sample sample1 \
    --fq1 results/04_bwa/sample1/sample1_species_1.fq \
    --fq2 results/04_bwa/sample1/sample1_species_2.fq \
    --reference /path/to/host_residual.fa \
    --outdir results/05_host_filter/sample1 \
    --db /path/to/kraken2_db

Optional:
  --bwa, --samtools, --seqkit, --kraken2, --rscript, --threads
EOF
}

SAMPLE=""
FQ1=""
FQ2=""
REFERENCE=""
OUTDIR=""
DB=""
BWA="bwa"
SAMTOOLS="samtools"
SEQKIT="seqkit"
KRAKEN2="kraken2"
RSCRIPT="Rscript"
THREADS=8
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --fq1) FQ1="$2"; shift 2 ;;
    --fq2) FQ2="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --db) DB="$2"; shift 2 ;;
    --bwa) BWA="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    --seqkit) SEQKIT="$2"; shift 2 ;;
    --kraken2) KRAKEN2="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLE" && -n "$FQ1" && -n "$FQ2" && -n "$REFERENCE" && -n "$OUTDIR" && -n "$DB" ]] || {
  usage >&2
  exit 1
}
[[ -s "$FQ1" && -s "$FQ2" ]] || { echo "Input FASTQ not found" >&2; exit 1; }
[[ -s "${REFERENCE}.bwt" ]] || {
  echo "BWA index not found for $REFERENCE; run: bwa index $REFERENCE" >&2
  exit 1
}
mkdir -p "$OUTDIR"

echo "[05] ${SAMPLE}: host-residual filtering"
"$RSCRIPT" "$SCRIPT_DIR/05_filter_host_residuals.R" \
  --sample "$SAMPLE" --fq1 "$FQ1" --fq2 "$FQ2" \
  --reference "$REFERENCE" --outdir "$OUTDIR" \
  --bwa "$BWA" --samtools "$SAMTOOLS" --seqkit "$SEQKIT" --threads "$THREADS"

echo "[05] ${SAMPLE}: refreshing read-to-taxid assignments"
bash "$SCRIPT_DIR/utils/refresh_taxid_assignments.sh" \
  --sample "$SAMPLE" \
  --r1 "$OUTDIR/${SAMPLE}_clean_1.fq" \
  --r2 "$OUTDIR/${SAMPLE}_clean_2.fq" \
  --outdir "$OUTDIR" --db "$DB" \
  --kraken2 "$KRAKEN2" --threads "$THREADS" \
  --rscript "$RSCRIPT" --script-dir "$SCRIPT_DIR"

echo "[05] done: $SAMPLE"
