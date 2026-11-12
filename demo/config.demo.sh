#!/usr/bin/env bash
# Demo configuration. Run from the repository root:
#   bash scripts/run_pipeline.sh --config demo/config.demo.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="$PROJECT_ROOT/demo/samples.tsv"
OUTDIR="$PROJECT_ROOT/demo/results"

RSCRIPT="Rscript"
PYTHON="python3"
KRAKEN2_BIN="kraken2"
BWA="bwa"
SAMTOOLS="samtools"
BAMTOFASTQ="$PROJECT_ROOT/demo/utils/bamtofastq_samtools.sh"
SEQKIT="seqkit"
SEQTK="seqtk"
TAXONKIT="taxonkit"

KRAKEN_DB="$PROJECT_ROOT/demo/reference/kraken_db"
BWA_LIBRARY="$PROJECT_ROOT/demo/reference/bwa_library"
TAXON_LINEAGE="$PROJECT_ROOT/demo/reference/taxonkit_lineage_matrix.tsv"
CURATED_SPECIES="$PROJECT_ROOT/demo/reference/curated_bacteria_species.tsv"
HOST_RESIDUAL_FASTA="$PROJECT_ROOT/demo/reference/host_residual.fa"
CELLLINE_Q90="$PROJECT_ROOT/demo/reference/cell_line_rpmm_q90.tsv"
VIRUS_NCBI_HOST="$PROJECT_ROOT/demo/reference/ncbi_virus_host_demo.tsv"
VIRUS_NCBI_ASSEMBLY_HOST=""
VIRUS_HVD="$PROJECT_ROOT/demo/reference/hvd_demo.csv"
VIRUS_HVIDB="$PROJECT_ROOT/demo/reference/hvidb_demo.txt"
VIRUS_MIN_SCORE=2
KEEP_OTHER_KINGDOMS=0
BLACK_GENUS="$PROJECT_ROOT/demo/reference/black_list_demo.txt"
MGNIFY_BLACK="$PROJECT_ROOT/demo/reference/mgnify_black_demo.csv"
CCLE_BLACKLIST="$PROJECT_ROOT/demo/reference/ccle_blacklist_demo.txt"

BC_LEN=16
UMI_LEN=10
THREADS=4
BWA_THREADS=2
TAXA_THREADS=2

MIN_SCORE=6
MIN_ALN=45
MIN_IDENTITY=0.95
MIN_KMER_FRACTION=0.5
MIN_CORR=0.1
DECONTAM_THRESHOLD=0.1
MGNIFY_MIN_STUDIES=4
CELLRANGER_FILTER=1
