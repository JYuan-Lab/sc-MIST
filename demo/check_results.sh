#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED="$ROOT/demo/expected"
RESULTS="$ROOT/demo/results"
FAIL=0

TMP="$(mktemp -d "${TMPDIR:-/tmp}/scmist_demo.XXXXXX")"
trap 'find "$TMP" -type f -delete 2>/dev/null; rmdir "$TMP" 2>/dev/null || true' EXIT

check() {
  local name="$1" expected="$2" observed="$3"
  if [[ ! -s "$observed" ]]; then
    echo "[FAIL] missing result: $observed" >&2
    FAIL=1
  elif diff -u "$expected" "$observed" >"$TMP/diff" 2>&1; then
    echo "[ok]   $name"
  else
    echo "[FAIL] $name differs from expected:" >&2
    sed -n '1,20p' "$TMP/diff" >&2
    FAIL=1
  fi
}

# 11_summary contains the two project-level tables.
check "taxon x sample matrix" \
  "$EXPECTED/taxon_umi_matrix.tsv" \
  "$RESULTS/11_summary/taxon_umi_matrix.tsv"
check "cell x taxon long table" \
  "$EXPECTED/cell_taxon_umi_long.tsv" \
  "$RESULTS/11_summary/cell_taxon_umi_long.tsv"

# Build the QC summary in the same order as the expected file.
{
  printf 'sample\tmicrobial_umi\ttaxa_detected\tbarcodes_with_microbes\n'
  for s in sample1 sample2 sample3; do
    tail -n +2 "$RESULTS/10_quantify/$s/$s.qc.tsv"
  done
} > "$TMP/qc_summary.tsv"

check "per-sample QC summary" "$EXPECTED/qc_summary.tsv" "$TMP/qc_summary.tsv"

# The decontam table must have been read by step 11 and must not reject a taxon.
if [[ ! -s "$RESULTS/09_decontam/decontam.tsv" ]]; then
  echo "[FAIL] missing decontam.tsv" >&2
  FAIL=1
elif ! awk -F'\t' 'NR==1 {for(j=1;j<=NF;j++) if($j=="p") pi=j; next}
                    pi && $pi!="NA" && $pi!="" && $pi<=0.1 {bad=1}
                    END {exit !(pi && !bad)}' "$RESULTS/09_decontam/decontam.tsv"; then
  echo "[FAIL] decontam.tsv has no usable p column or rejects every taxon" >&2
  FAIL=1
else
  echo "[ok]   decontam p-values"
fi

if [[ "$FAIL" -ne 0 ]]; then
  echo "Demo comparison failed." >&2
  exit 1
fi
echo "Demo comparison passed."
