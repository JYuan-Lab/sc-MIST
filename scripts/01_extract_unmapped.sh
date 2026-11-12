#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Extract human-unmapped read pairs from a Cell Ranger BAM and rebuild the
10x R1 read as barcode + UMI.

Required:
  --sample        Sample ID
  --bam           Cell Ranger possorted_genome_bam.bam
  --outdir        Output directory (one per sample)
  --bc-len        Cell barcode length (10x v2: 16, v3: 16)
  --umi-len       UMI length (10x v2: 10, v3: 12)

Optional:
  --threads       Threads for bamtofastq (default: 8)
  --samtools      samtools executable (default: samtools)
  --bamtofastq    bamtofastq executable (default: bamtofastq)
  --seqkit        seqkit executable (default: seqkit)
EOF
}

SAMPLE=""
BAM=""
OUTDIR=""
BC_LEN=""
UMI_LEN=""
THREADS=8
SAMTOOLS="samtools"
BAMTOFASTQ="bamtofastq"
SEQKIT="seqkit"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bc-len) BC_LEN="$2"; shift 2 ;;
    --umi-len) UMI_LEN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    --bamtofastq) BAMTOFASTQ="$2"; shift 2 ;;
    --seqkit) SEQKIT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLE" && -n "$BAM" && -n "$OUTDIR" && -n "$BC_LEN" && -n "$UMI_LEN" ]] || {
  usage >&2
  exit 1
}
[[ -s "$BAM" ]] || { echo "BAM not found: $BAM" >&2; exit 1; }
command -v "$SAMTOOLS" >/dev/null || { echo "samtools not found: $SAMTOOLS" >&2; exit 1; }
command -v "$BAMTOFASTQ" >/dev/null || { echo "bamtofastq not found: $BAMTOFASTQ" >&2; exit 1; }
command -v "$SEQKIT" >/dev/null || { echo "seqkit not found: $SEQKIT" >&2; exit 1; }

BCUMI_LEN=$((BC_LEN + UMI_LEN))
mkdir -p "$OUTDIR"
UNMAPPED_BAM="$OUTDIR/${SAMPLE}.unmapped.bam"
TMPDIR="$OUTDIR/${SAMPLE}.bamtofastq"

echo "[01] ${SAMPLE}: extracting reads with flag 4 (unmapped)"
"$SAMTOOLS" view -b -f 4 "$BAM" > "$UNMAPPED_BAM"

mkdir -p "$TMPDIR"
echo "[01] ${SAMPLE}: converting unmapped BAM back to FASTQ"
"$BAMTOFASTQ" --nthreads "$THREADS" "$UNMAPPED_BAM" "$TMPDIR"

R1_FILES=()
while IFS= read -r f; do R1_FILES+=("$f"); done < <(
  find "$TMPDIR" -type f \( -name '*R1*.fastq.gz' -o -name '*R1*.fq.gz' \) | sort
)
R2_FILES=()
while IFS= read -r f; do R2_FILES+=("$f"); done < <(
  find "$TMPDIR" -type f \( -name '*R2*.fastq.gz' -o -name '*R2*.fq.gz' \) | sort
)

[[ ${#R1_FILES[@]} -gt 0 ]] || { echo "bamtofastq produced no R1 FASTQ under $TMPDIR" >&2; exit 1; }
[[ ${#R2_FILES[@]} -gt 0 ]] || { echo "bamtofastq produced no R2 FASTQ under $TMPDIR" >&2; exit 1; }
[[ ${#R1_FILES[@]} -eq ${#R2_FILES[@]} ]] || {
  echo "Mismatched R1/R2 FASTQ counts under $TMPDIR" >&2
  exit 1
}

R1="$OUTDIR/${SAMPLE}_R1.fastq.gz"
R2="$OUTDIR/${SAMPLE}_R2.fastq.gz"
R1_BCUMI="$OUTDIR/${SAMPLE}_R1.bcumi.fastq.gz"

cat "${R1_FILES[@]}" > "$R1"
cat "${R2_FILES[@]}" > "$R2"

echo "[01] ${SAMPLE}: keeping the first ${BCUMI_LEN} bp of R1 (barcode ${BC_LEN} + UMI ${UMI_LEN})"
"$SEQKIT" subseq -r "1:${BCUMI_LEN}" "$R1" | gzip -c > "$R1_BCUMI"

find "$TMPDIR" -type f -delete
find "$TMPDIR" -depth -type d -empty -delete

cat <<EOF
[01] done: ${SAMPLE}
  ${UNMAPPED_BAM}
  ${R1}
  ${R2}
  ${R1_BCUMI}
EOF
