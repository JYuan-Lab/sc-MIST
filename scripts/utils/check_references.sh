#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Check that the reference files required by sc-MIST are present and have the
expected basic structure. This does not validate biological content.

Usage:
  bash scripts/utils/check_references.sh --config config/config.sh
  bash scripts/utils/check_references.sh --config demo/config.demo.sh
USAGE
}

CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done
[[ -s "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; usage >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

FAIL=0
ok()   { echo "[ok]   $*"; }
bad()  { echo "[FAIL] $*" >&2; FAIL=1; }

check_file() {
  local label="$1" path="$2"
  if [[ -s "$path" ]]; then ok "$label: $path"; else bad "$label missing or empty: $path"; fi
}

check_header() {
  local label="$1" path="$2" expected="$3"
  [[ -s "$path" ]] || { bad "$label missing: $path"; return; }
  local header
  header="$(head -n 1 "$path")"
  if [[ "$header" == *"$expected"* ]]; then
    ok "$label header"
  else
    bad "$label header does not contain '$expected' (found: $header)"
  fi
}

echo "Checking sc-MIST references"
check_file "Kraken2 hash"   "${KRAKEN_DB}/hash.k2d"
check_file "Kraken2 opts"   "${KRAKEN_DB}/opts.k2d"
check_file "Kraken2 taxo"   "${KRAKEN_DB}/taxo.k2d"
check_file "HLA/rRNA FASTA" "$HOST_RESIDUAL_FASTA"
check_file "HLA/rRNA BWA index" "${HOST_RESIDUAL_FASTA}.bwt"
check_file "lineage matrix" "$TAXON_LINEAGE"
check_file "curated species" "$CURATED_SPECIES"
check_file "cell-line q90"  "$CELLLINE_Q90"
check_file "blacklist"      "$BLACK_GENUS"
check_file "MGnify blacklist" "$MGNIFY_BLACK"
check_file "CCLE blacklist" "$CCLE_BLACKLIST"
check_file "HVD virus-organ table" "${VIRUS_HVD:-}"
check_file "HVIDB interaction table" "${VIRUS_HVIDB:-}"
if [[ -n "${VIRUS_NCBI_HOST:-}" ]]; then
  check_file "NCBI Virus host metadata" "$VIRUS_NCBI_HOST"
fi

check_header "lineage matrix" "$TAXON_LINEAGE" $'taxid\tkindom\tphylum\tgenus\tspecies\tg_taxid\ts_taxid'
check_header "curated species" "$CURATED_SPECIES" $'species\tspecies_taxid'
if head -n 1 "$CELLLINE_Q90" 2>/dev/null | grep -q 'threshold_0.9_exclude0'; then
  check_header "cell-line q90" "$CELLLINE_Q90" 'taxid'
else
  check_header "cell-line q90" "$CELLLINE_Q90" $'taxid\tq90'
fi

check_header "HVD table" "${VIRUS_HVD:-}" 'VirusName'
check_header "HVIDB table" "${VIRUS_HVIDB:-}" 'Organism_Interactor_virus'
if [[ -n "${VIRUS_NCBI_HOST:-}" ]]; then
  check_header "NCBI Virus host metadata: Organism_Name" "$VIRUS_NCBI_HOST" 'Organism_Name'
  check_header "NCBI Virus host metadata: Species" "$VIRUS_NCBI_HOST" 'Species'
  check_header "NCBI Virus host metadata: Genus" "$VIRUS_NCBI_HOST" 'Genus'
  check_header "NCBI Virus host metadata: Host" "$VIRUS_NCBI_HOST" 'Host'
fi

if [[ -d "$BWA_LIBRARY" ]]; then
  n_bwt=$(find "$BWA_LIBRARY" -maxdepth 1 -type f -name '*.fa.bwt' | wc -l | tr -d ' ')
  if [[ "$n_bwt" -gt 0 ]]; then
    ok "BWA library: $n_bwt indexed FASTA file(s) in $BWA_LIBRARY"
  else
    bad "BWA library has no *.fa.bwt files: $BWA_LIBRARY"
  fi
else
  bad "BWA library directory missing: $BWA_LIBRARY"
fi

if [[ "$FAIL" -ne 0 ]]; then
  echo "Reference check failed." >&2
  exit 1
fi
echo "Reference check passed."
