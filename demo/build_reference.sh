#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$ROOT/demo/reference"
THREADS="${THREADS:-2}"

command -v kraken2-build >/dev/null || { echo "kraken2-build not found in PATH" >&2; exit 1; }
command -v bwa >/dev/null || { echo "bwa not found in PATH" >&2; exit 1; }

DB="$REF/kraken_db"
if [[ ! -s "$DB/hash.k2d" ]]; then
  echo "[demo] building toy Kraken2 database"
  mkdir -p "$DB/taxonomy"
  cp "$REF/kraken_taxonomy/nodes.dmp" "$REF/kraken_taxonomy/names.dmp" "$DB/taxonomy/"
  mkdir -p "$DB/library"
  find "$DB/library" -type f -delete
  find "$DB/library" -depth -type d -empty -delete
  mkdir -p "$DB/library"
  for f in "$REF"/taxon_*.fa; do
    kraken2-build --add-to-library "$f" --db "$DB"
  done
  kraken2-build --build --db "$DB" --threads "$THREADS" \
    --kmer-len 35 --minimizer-len 31 --no-masking
else
  echo "[demo] toy Kraken2 database already exists"
fi

BWA_LIB="$REF/bwa_library"
mkdir -p "$BWA_LIB"

make_ref() {
  local taxid="$1"
  shift
  local out="$BWA_LIB/${taxid}.fa"
  : > "$out"
  for src in "$@"; do
    echo ">taxid_${taxid}_${src##*/}" >> "$out"
    tail -n +2 "$src" >> "$out"
  done
  if [[ ! -s "$out.bwt" ]]; then
    bwa index "$out"
  fi
}

make_ref 1001 "$REF/taxon_1001.fa"
make_ref 1002 "$REF/taxon_1002.fa"
make_ref 1000 "$REF/taxon_1001.fa" "$REF/taxon_1002.fa"
make_ref 2001 "$REF/taxon_2001.fa"
make_ref 2000 "$REF/taxon_2001.fa"
make_ref 3001 "$REF/taxon_3001.fa"
make_ref 3000 "$REF/taxon_3001.fa"

if [[ ! -s "$REF/host_residual.fa.bwt" ]]; then
  bwa index "$REF/host_residual.fa"
fi

echo "[demo] reference indices are ready under $REF"
