#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build the taxon-specific BWA reference library used by step 04.

Input species FASTAs must be named <taxid>.fa (or .fasta/.fna). The lineage
matrix must contain the columns: taxid, g_taxid, s_taxid.

Usage:
  bash scripts/utils/build_bwa_library.sh \
    --species-dir /path/to/species_fasta \
    --lineage /path/to/taxonkit_lineage_matrix.tsv \
    --outdir /path/to/bwa_library
EOF
}

SPECIES_DIR=""
LINEAGE=""
OUTDIR=""
BWA="bwa"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species-dir) SPECIES_DIR="$2"; shift 2 ;;
    --lineage) LINEAGE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bwa) BWA="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "$SPECIES_DIR" && -s "$LINEAGE" && -n "$OUTDIR" ]] || { usage >&2; exit 1; }
command -v "$BWA" >/dev/null || { echo "bwa not found: $BWA" >&2; exit 1; }
mkdir -p "$OUTDIR"

find_species_fa() {
  local taxid="$1"
  for ext in fa fasta fna; do
    if [[ -s "$SPECIES_DIR/${taxid}.${ext}" ]]; then
      echo "$SPECIES_DIR/${taxid}.${ext}"
      return 0
    fi
  done
  return 1
}

for f in "$SPECIES_DIR"/*.fa "$SPECIES_DIR"/*.fasta "$SPECIES_DIR"/*.fna; do
  [[ -s "$f" ]] || continue
  taxid="$(basename "$f")"
  taxid="${taxid%%.*}"
  out="$OUTDIR/${taxid}.fa"
  [[ -s "$out" ]] || cp "$f" "$out"
  [[ -s "$out.bwt" ]] || "$BWA" index "$out"
done

PAIRS="$(mktemp)"
trap 'find "$PAIRS" -type f -delete' EXIT
awk -F'\t' 'NR > 1 && $6 != "NA" && $7 != "NA" && $6 != "" && $7 != "" {print $6 "\t" $7}' \
  "$LINEAGE" | sort -u > "$PAIRS"

cut -f1 "$PAIRS" | sort -u | while read -r genus; do
  out="$OUTDIR/${genus}.fa"
  tmp="$OUTDIR/.${genus}.fa.tmp"
  : > "$tmp"
  awk -F'\t' -v g="$genus" '$1 == g {print $2}' "$PAIRS" | while read -r sp; do
    src="$(find_species_fa "$sp" || true)"
    if [[ -n "$src" ]]; then
      cat "$src" >> "$tmp"
    else
      echo "WARNING: no species FASTA for taxid $sp (genus $genus)" >&2
    fi
  done
  if [[ -s "$tmp" ]]; then
    mv "$tmp" "$out"
    [[ -s "$out.bwt" ]] || "$BWA" index "$out"
  else
    find "$tmp" -type f -delete
  fi
done

echo "[bwa-library] done: $OUTDIR"
