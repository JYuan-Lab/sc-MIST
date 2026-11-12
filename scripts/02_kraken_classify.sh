#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Classify paired reads with Kraken2/Kraken2Uniq and create a standard Kraken
report plus an MPA-style report.

Required:
  --sample          Sample ID
  --r1              R1 FASTQ (barcode + UMI, optionally gzipped)
  --r2              R2 FASTQ (optionally gzipped)
  --outdir          Output directory
  --db              Kraken2 database directory

Optional:
  --kraken2         kraken2 or Kraken2Uniq executable (default: kraken2)
  --threads         Classification threads (default: 8)
  --kreport2mpa     kreport2mpa.py path (default: scripts/utils/kreport2mpa.py)
EOF
}

SAMPLE=""
R1=""
R2=""
OUTDIR=""
DB=""
KRAKEN2="kraken2"
THREADS=8
KREPORT2MPA="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/utils/kreport2mpa.py"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --db) DB="$2"; shift 2 ;;
    --kraken2) KRAKEN2="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --kreport2mpa) KREPORT2MPA="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLE" && -n "$R1" && -n "$R2" && -n "$OUTDIR" && -n "$DB" ]] || {
  usage >&2
  exit 1
}
[[ -s "$R1" && -s "$R2" ]] || { echo "R1/R2 not found" >&2; exit 1; }
command -v "$KRAKEN2" >/dev/null || { echo "Kraken2 not found: $KRAKEN2" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 1; }
[[ -s "$KREPORT2MPA" ]] || { echo "kreport2mpa.py not found: $KREPORT2MPA" >&2; exit 1; }

mkdir -p "$OUTDIR"
OUT="$OUTDIR/${SAMPLE}.kraken.output.txt"
REPORT="$OUTDIR/${SAMPLE}.kraken.report.txt"
REPORT_STD="$OUTDIR/${SAMPLE}.kraken.report.std.txt"
MPA="$OUTDIR/${SAMPLE}.kraken.report.mpa.txt"

echo "[02] ${SAMPLE}: Kraken2 classification"
"$KRAKEN2" \
  --db "$DB" \
  --threads "$THREADS" \
  --paired \
  --use-names \
  --report-minimizer-data \
  --classified-out "${OUTDIR}/${SAMPLE}#.fq" \
  --output "$OUT" \
  --report "$REPORT" \
  "$R1" "$R2"

# Keep the columns used by the original downstream scripts, without the
# minimizer columns added by --report-minimizer-data.
cut -f1-3,6-8 "$REPORT" > "$REPORT_STD"
python3 "$KREPORT2MPA" -r "$REPORT_STD" -o "$MPA" --intermediate-ranks

cat <<EOF
[02] done: ${SAMPLE}
  ${OUT}
  ${REPORT}
  ${REPORT_STD}
  ${MPA}
  ${OUTDIR}/${SAMPLE}_1.fq
  ${OUTDIR}/${SAMPLE}_2.fq
EOF
