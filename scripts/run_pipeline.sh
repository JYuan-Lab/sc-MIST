#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$ROOT/config/config.sh"
SAMPLES_TSV=""
OUTDIR=""
CLI_SAMPLES=""
CLI_OUTDIR=""
START=1
END=11
SKIP_DECONTAM=0

usage() {
  cat <<'EOF'
Run the sc-MIST pipeline for all samples listed in a samples TSV.

Usage:
  bash scripts/run_pipeline.sh [options]

Options:
  --config FILE          Config file (default: config/config.sh)
  --samples FILE         samples.tsv (default: from config)
  --outdir DIR           Output directory (default: from config)
  --start N              First step to run (1-11, default 1)
  --end N                Last step to run (1-11, default 11)
  --skip-decontam        Write NA decontam p-values (for demos/tests)
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --samples) CLI_SAMPLES="$2"; shift 2 ;;
    --outdir) CLI_OUTDIR="$2"; shift 2 ;;
    --start) START="$2"; shift 2 ;;
    --end) END="$2"; shift 2 ;;
    --skip-decontam) SKIP_DECONTAM=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -s "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"
[[ -n "$CLI_SAMPLES" ]] && SAMPLES_TSV="$CLI_SAMPLES"
[[ -n "$CLI_OUTDIR" ]] && OUTDIR="$CLI_OUTDIR"
[[ -n "$SAMPLES_TSV" ]] || SAMPLES_TSV="$ROOT/config/samples.tsv"
[[ -n "$OUTDIR" ]] || OUTDIR="$ROOT/results"
[[ -s "$SAMPLES_TSV" ]] || { echo "samples.tsv not found: $SAMPLES_TSV" >&2; exit 1; }

DIRS=(01_unmapped 02_kraken 03_microbial 04_bwa 05_host_filter \
      06_denoise 07_sample_denoise 08_cellline 09_decontam 10_quantify 11_summary)
mkdir -p "$OUTDIR/logs"
for d in "${DIRS[@]}"; do mkdir -p "$OUTDIR/$d"; done

SAMPLE_IDS=()
while IFS=$'\t' read -r SID _REST; do
  [[ -n "$SID" && "$SID" != "sample_id" ]] && SAMPLE_IDS+=("$SID")
done < <(tail -n +2 "$SAMPLES_TSV")
if [[ ${#SAMPLE_IDS[@]} -lt 3 ]]; then
  echo "WARNING: fewer than 3 samples; sample-level correlation tests require >=3 samples." >&2
fi
SAMPLES_COMMA=$(IFS=,; echo "${SAMPLE_IDS[*]}")

run_step() {
  local step="$1"
  if (( step >= START && step <= END )); then return 0; else return 1; fi
}

for S in "${SAMPLE_IDS[@]}"; do
  CR_OUT=$(awk -F'\t' -v s="$S" '$1==s {print $2; exit}' "$SAMPLES_TSV")
  [[ -n "$CR_OUT" ]] || { echo "No Cell Ranger outs path for $S" >&2; exit 1; }
  [[ "$CR_OUT" = /* ]] || CR_OUT="$ROOT/$CR_OUT"
  for d in 01_unmapped 02_kraken 03_microbial 04_bwa 05_host_filter 06_denoise 08_cellline 10_quantify; do
    mkdir -p "$OUTDIR/$d/$S"
  done

  if run_step 1; then
    bash "$ROOT/scripts/01_extract_unmapped.sh" \
      --sample "$S" --bam "$CR_OUT/possorted_genome_bam.bam" \
      --outdir "$OUTDIR/01_unmapped/$S" \
      --bc-len "$BC_LEN" --umi-len "$UMI_LEN" --threads "$THREADS" \
      --samtools "$SAMTOOLS" --bamtofastq "$BAMTOFASTQ" --seqkit "$SEQKIT" \
      2>&1 | tee "$OUTDIR/logs/${S}.01.log"
  fi

  if run_step 2; then
    bash "$ROOT/scripts/02_kraken_classify.sh" \
      --sample "$S" \
      --r1 "$OUTDIR/01_unmapped/$S/${S}_R1.bcumi.fastq.gz" \
      --r2 "$OUTDIR/01_unmapped/$S/${S}_R2.fastq.gz" \
      --outdir "$OUTDIR/02_kraken/$S" --db "$KRAKEN_DB" \
      --kraken2 "$KRAKEN2_BIN" --threads "$THREADS" \
      2>&1 | tee "$OUTDIR/logs/${S}.02.log"
  fi

  if run_step 3; then
    "$RSCRIPT" "$ROOT/scripts/03_extract_microbial_reads.R" \
      --sample "$S" \
      --kraken-report "$OUTDIR/02_kraken/$S/${S}.kraken.report.txt" \
      --mpa-report "$OUTDIR/02_kraken/$S/${S}.kraken.report.mpa.txt" \
      --kraken-output "$OUTDIR/02_kraken/$S/${S}.kraken.output.txt" \
      --fq1 "$OUTDIR/02_kraken/$S/${S}_1.fq" \
      --fq2 "$OUTDIR/02_kraken/$S/${S}_2.fq" \
      --outdir "$OUTDIR/03_microbial/$S" \
      2>&1 | tee "$OUTDIR/logs/${S}.03.log"
  fi

  if run_step 4; then
    "$RSCRIPT" "$ROOT/scripts/04_validate_taxa_bwa.R" \
      --sample "$S" \
      --microbiome-output "$OUTDIR/03_microbial/$S/${S}.microbiome.output.txt" \
      --fa1 "$OUTDIR/03_microbial/$S/${S}_1.fa" \
      --fa2 "$OUTDIR/03_microbial/$S/${S}_2.fa" \
      --lineage "$TAXON_LINEAGE" --curated-species "$CURATED_SPECIES" \
      --bwa-library "$BWA_LIBRARY" --outdir "$OUTDIR/04_bwa/$S" \
      --bwa "$BWA" --samtools "$SAMTOOLS" --seqkit "$SEQKIT" --seqtk "$SEQTK" \
      --bwa-filter "$ROOT/scripts/utils/bwa_filter.awk" \
      --min-score "$MIN_SCORE" --min-aln "$MIN_ALN" --min-identity "$MIN_IDENTITY" \
      --threads "$TAXA_THREADS" --bwa-threads "$BWA_THREADS" \
      2>&1 | tee "$OUTDIR/logs/${S}.04.log"
  fi

  if run_step 5; then
    bash "$ROOT/scripts/05_filter_host_and_refresh.sh" \
      --sample "$S" \
      --fq1 "$OUTDIR/04_bwa/$S/${S}_species_1.fq" \
      --fq2 "$OUTDIR/04_bwa/$S/${S}_species_2.fq" \
      --reference "$HOST_RESIDUAL_FASTA" --outdir "$OUTDIR/05_host_filter/$S" \
      --db "$KRAKEN_DB" --bwa "$BWA" --samtools "$SAMTOOLS" --seqkit "$SEQKIT" \
      --kraken2 "$KRAKEN2_BIN" --rscript "$RSCRIPT" \
      --threads "$THREADS" --script-dir "$ROOT/scripts" \
      2>&1 | tee "$OUTDIR/logs/${S}.05.log"
  fi

  if run_step 6; then
    "$RSCRIPT" "$ROOT/scripts/06_barcode_denoising.R" \
      --sample "$S" \
      --fa1 "$OUTDIR/05_host_filter/$S/microbial/${S}_1.fa" \
      --fa2 "$OUTDIR/05_host_filter/$S/microbial/${S}_2.fa" \
      --kraken-report "$OUTDIR/05_host_filter/$S/${S}.kraken.report.txt" \
      --mpa-report "$OUTDIR/05_host_filter/$S/${S}.kraken.report.mpa.txt" \
      --microbiome-output "$OUTDIR/05_host_filter/$S/microbial/${S}.microbiome.output.txt" \
      --first-kraken-report "$OUTDIR/02_kraken/$S/${S}.kraken.report.txt" \
      --outdir "$OUTDIR/06_denoise/$S" \
      --cb-len "$BC_LEN" --umi-len "$UMI_LEN" --min-frac "$MIN_KMER_FRACTION" \
      2>&1 | tee "$OUTDIR/logs/${S}.06.log"
  fi

done

if run_step 7; then
  "$RSCRIPT" "$ROOT/scripts/07_sample_denoising.R" \
    --samples "$SAMPLES_COMMA" --rpmm-dir "$OUTDIR/06_denoise" \
    --outdir "$OUTDIR/07_sample_denoise" \
    2>&1 | tee "$OUTDIR/logs/07_sample_denoising.log"
fi

for S in "${SAMPLE_IDS[@]}"; do
  if run_step 8; then
    "$RSCRIPT" "$ROOT/scripts/08_filter_cellline.R" \
      --sample "$S" \
      --barcode-hits "$OUTDIR/06_denoise/$S/${S}.barcode.kmer.hits.tsv" \
      --rpmm "$OUTDIR/06_denoise/$S/${S}.kraken.report.rpmm.tsv" \
      --sample-denoising "$OUTDIR/07_sample_denoise/sample_denoising.tsv" \
      --kraken-report "$OUTDIR/05_host_filter/$S/${S}.kraken.report.txt" \
      --cellline-q90 "$CELLLINE_Q90" --outdir "$OUTDIR/08_cellline/$S" \
      2>&1 | tee "$OUTDIR/logs/${S}.08.log"
  fi
done

if run_step 9; then
  EXTRA=()
  if [[ "$SKIP_DECONTAM" -eq 1 ]]; then EXTRA+=(--skip-decontam); fi
  "$RSCRIPT" "$ROOT/scripts/09_decontam.R" \
    --samples "$SAMPLES_COMMA" --hits-dir "$OUTDIR/08_cellline" \
    --kraken-dir "$OUTDIR/02_kraken" --outdir "$OUTDIR/09_decontam" \
    --threshold "$DECONTAM_THRESHOLD" ${EXTRA[@]+"${EXTRA[@]}"} \
    2>&1 | tee "$OUTDIR/logs/09_decontam.log"
fi

for S in "${SAMPLE_IDS[@]}"; do
  CR_OUT=$(awk -F'\t' -v s="$S" '$1==s {print $2; exit}' "$SAMPLES_TSV")
  [[ -n "$CR_OUT" ]] || { echo "No Cell Ranger outs path for $S" >&2; exit 1; }
  [[ "$CR_OUT" = /* ]] || CR_OUT="$ROOT/$CR_OUT"

  if run_step 10; then
    EXTRA=()
    if [[ "$CELLRANGER_FILTER" -eq 1 ]]; then
      EXTRA+=(--cellranger-outs "$CR_OUT")
    else
      EXTRA+=(--no-cellranger-filter)
    fi
    "$RSCRIPT" "$ROOT/scripts/10_quantify_filter.R" \
      --sample "$S" \
      --fa1 "$OUTDIR/05_host_filter/$S/microbial/${S}_1.fa" \
      --fa2 "$OUTDIR/05_host_filter/$S/microbial/${S}_2.fa" \
      --taxa "$OUTDIR/08_cellline/$S/${S}.cell_line_quantile_hits_taxa.tsv" \
      --kraken-report "$OUTDIR/05_host_filter/$S/${S}.kraken.report.txt" \
      --mpa-report "$OUTDIR/05_host_filter/$S/${S}.kraken.report.mpa.txt" \
      --hits "$OUTDIR/08_cellline/$S/${S}.cell_line_quantile_hits.tsv" \
      --decontam "$OUTDIR/09_decontam/decontam.tsv" \
      --bwa-black "$OUTDIR/04_bwa/$S/${S}_bwa_black_species_id.tsv" \
      --black-genus "$BLACK_GENUS" --mgnify "$MGNIFY_BLACK" \
      --ccle-blacklist "$CCLE_BLACKLIST" \
      --outdir "$OUTDIR/10_quantify/$S" \
      --cb-len "$BC_LEN" --umi-len "$UMI_LEN" \
      --min-corr "$MIN_CORR" --decontam-threshold "$DECONTAM_THRESHOLD" \
      --mgnify-min-studies "$MGNIFY_MIN_STUDIES" \
      ${EXTRA[@]+"${EXTRA[@]}"} \
      2>&1 | tee "$OUTDIR/logs/${S}.10.log"
  fi
done

if run_step 11; then
  VIRUS_ARGS=(--lineage "$TAXON_LINEAGE" --virus-min-score "$VIRUS_MIN_SCORE")
  if [[ -n "${VIRUS_NCBI_HOST:-}" ]]; then VIRUS_ARGS+=(--ncbi-virus-host "$VIRUS_NCBI_HOST"); fi
  if [[ -n "${VIRUS_NCBI_ASSEMBLY_HOST:-}" ]]; then VIRUS_ARGS+=(--ncbi-virus-assembly-host "$VIRUS_NCBI_ASSEMBLY_HOST"); fi
  if [[ -n "${VIRUS_HVD:-}" ]]; then VIRUS_ARGS+=(--hvd "$VIRUS_HVD"); fi
  if [[ -n "${VIRUS_HVIDB:-}" ]]; then VIRUS_ARGS+=(--hvidb "$VIRUS_HVIDB"); fi
  if [[ "${KEEP_OTHER_KINGDOMS:-0}" -eq 1 ]]; then VIRUS_ARGS+=(--keep-other-kingdoms); fi
  "$RSCRIPT" "$ROOT/scripts/11_summarize.R" \
    --samples "$SAMPLES_COMMA" \
    --abundance-dir "$OUTDIR/10_quantify" --counts-dir "$OUTDIR/10_quantify" \
    --outdir "$OUTDIR/11_summary" \
    ${VIRUS_ARGS[@]+"${VIRUS_ARGS[@]}"} \
    2>&1 | tee "$OUTDIR/logs/11_summary.log"
fi

echo "sc-MIST pipeline finished. Results: $OUTDIR"
