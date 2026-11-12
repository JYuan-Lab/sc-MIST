# Keep primary alignments with alignment length > min_len and identity > min_id.
# Input is SAM from `samtools view -F 0x900`.
BEGIN { FS = "\t" }
{
  cigar = $6
  nm = 0
  aln = 0
  for (i = 12; i <= NF; i++) {
    if ($i ~ /^NM:i:/) {
      split($i, a, ":")
      nm = a[3]
    }
  }
  tmp = cigar
  while (match(tmp, /[0-9]+[M=X]/)) {
    aln += substr(tmp, RSTART, RLENGTH - 1)
    tmp = substr(tmp, RSTART + RLENGTH)
  }
  if (aln > min_len && (1 - nm / aln) > min_id) {
    print $1
  }
}
