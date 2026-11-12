# Minimal FASTA reader used by steps 07 and 11.
# Returns a named character vector: names are headers without ">", values are
# sequences with line breaks removed. Keeping this in base R avoids a heavy
# Bioconductor dependency for two small parsing calls.
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) return(setNames(character(0), character(0)))
  is_header <- startsWith(lines, ">")
  headers <- sub("^>", "", lines[is_header])
  if (length(headers) == 0) return(setNames(character(0), character(0)))
  grp <- cumsum(is_header)
  seqs <- tapply(
    lines[!is_header],
    factor(grp[!is_header], levels = seq_along(headers)),
    paste0, collapse = ""
  )
  seqs[is.na(seqs)] <- ""
  out <- as.vector(seqs)
  names(out) <- headers
  out
}
