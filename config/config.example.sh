#!/usr/bin/env bash
# Copy this file to config/config.sh and edit the paths for your system.

# Project layout
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="$PROJECT_ROOT/config/samples.tsv"
OUTDIR="$PROJECT_ROOT/results"

# Executables
RSCRIPT="Rscript"
PYTHON="python3"
KRAKEN2_BIN="kraken2"
BWA="bwa"
SAMTOOLS="samtools"
BAMTOFASTQ="bamtofastq"
SEQKIT="seqkit"
SEQTK="seqtk"
TAXONKIT="taxonkit"

# Reference resources (see README "Reference preparation")
KRAKEN_DB="/path/to/kraken2_db"
BWA_LIBRARY="/path/to/bwa_library"
TAXON_LINEAGE="/path/to/taxonkit_lineage_matrix.tsv"
CURATED_SPECIES="$PROJECT_ROOT/resources/curated_bacteria_species.tsv"
HOST_RESIDUAL_FASTA="/path/to/host_residual.fa"
CELLLINE_Q90="$PROJECT_ROOT/resources/CCLE_quantile_mhc_rrna.txt"
BLACK_GENUS="$PROJECT_ROOT/resources/black_list.txt"
MGNIFY_BLACK="$PROJECT_ROOT/resources/mgnify_black.txt"
CCLE_BLACKLIST="$PROJECT_ROOT/resources/ccle_blacklist_taxids.txt"

# Viral host-evidence sources used by step 11.
# The NCBI Virus host metadata export is ~1.3 GB and is not bundled.
VIRUS_NCBI_HOST="/path/to/ncbi_virus_nucleotide_host_meta_1.tsv"
VIRUS_NCBI_ASSEMBLY_HOST=""   # optional second NCBI export
VIRUS_HVD="$PROJECT_ROOT/resources/humanvirusbase_Virus_Organ_pairs.csv"
VIRUS_HVIDB="$PROJECT_ROOT/resources/HVIDB_PPIs.txt"
VIRUS_MIN_SCORE=2             # keep viruses supported by >=2 of NCBI/HVD/HVIDB
KEEP_OTHER_KINGDOMS=0         # 1 = keep fungi/other kingdoms instead of bacteria-only

# 10x chemistry
BC_LEN=16
UMI_LEN=10

# Resource usage
THREADS=8
BWA_THREADS=4
TAXA_THREADS=2

# Filtering thresholds
MIN_SCORE=6
MIN_ALN=45
MIN_IDENTITY=0.95
MIN_KMER_FRACTION=0.5
MIN_CORR=0.1
DECONTAM_THRESHOLD=0.1
MGNIFY_MIN_STUDIES=4

# Set to 1 to restrict the final cell-level table to Cell Ranger barcodes.
CELLRANGER_FILTER=1
