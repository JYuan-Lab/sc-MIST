#!/usr/bin/env Rscript

# Build the seven-column taxonomy table used by step 04.
#
# The table is derived from an NCBI taxdump with taxonkit. It can be built
# either from a list of taxids or from the FASTA files in a BWA library.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--taxids", type = "character", default = NULL,
              help = "Text file with one taxid per line"),
  make_option("--bwa-library", dest = "bwa_library", type = "character", default = NULL,
              help = "Alternatively, derive taxids from <taxid>.fa files in this directory"),
  make_option("--taxdump", type = "character", help = "Directory containing nodes.dmp and names.dmp"),
  make_option("--out", type = "character", help = "Output TSV"),
  make_option("--taxonkit", type = "character", default = "taxonkit",
              help = "taxonkit executable [default %default]"),
  make_option("--lineage-raw", dest = "lineage_raw", type = "character", default = NULL,
              help = "Optional path to keep the raw taxonkit lineage output")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$taxids) == is.null(opt$bwa_library)) {
  stop("Provide exactly one of --taxids or --bwa-library")
}
if (is.null(opt$taxdump) || is.null(opt$out)) {
  stop("--taxdump and --out are required")
}
if (!dir.exists(opt$taxdump)) stop("taxdump directory not found: ", opt$taxdump)
if (nzchar(Sys.which(opt$taxonkit)) == FALSE) stop("taxonkit not found: ", opt$taxonkit)

if (!is.null(opt$bwa_library)) {
  if (!dir.exists(opt$bwa_library)) stop("BWA library not found: ", opt$bwa_library)
  fa <- list.files(opt$bwa_library, pattern = "\\.fa$", full.names = FALSE)
  taxids <- sub("\\.fa$", "", fa)
  taxids <- taxids[grepl("^[0-9]+$", taxids)]
} else {
  if (!file.exists(opt$taxids)) stop("taxid file not found: ", opt$taxids)
  taxids <- trimws(readLines(opt$taxids, warn = FALSE))
  taxids <- taxids[nzchar(taxids) & grepl("^[0-9]+$", taxids)]
}
taxids <- sort(unique(taxids))
if (length(taxids) == 0) stop("No numeric taxids found")

tmp_taxids <- tempfile(fileext = ".txt")
on.exit(unlink(tmp_taxids), add = TRUE)
writeLines(taxids, tmp_taxids)

raw <- if (is.null(opt$lineage_raw)) tempfile(fileext = ".tsv") else opt$lineage_raw
if (!is.null(opt$lineage_raw)) {
  dir.create(dirname(opt$lineage_raw), recursive = TRUE, showWarnings = FALSE)
}
status <- system2(
  opt$taxonkit,
  c("lineage", "--data-dir", opt$taxdump, "-t", "-R", tmp_taxids),
  stdout = raw, stderr = ""
)
if (status != 0) stop("taxonkit lineage failed with status ", status)

lin <- read.delim(raw, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE)
if (ncol(lin) < 4) stop("Unexpected taxonkit lineage output: need at least 4 columns")

rank_aliases <- list(
  kingdom = c("superkingdom", "kingdom", "domain"),
  phylum = c("phylum"),
  genus = c("genus"),
  species = c("species")
)
prefix <- c(kingdom = "k__", phylum = "p__", genus = "g__", species = "s__")

pick_rank <- function(names, taxids, ranks, level) {
  wanted <- rank_aliases[[level]]
  idx <- which(ranks %in% wanted)
  if (length(idx) == 0) return(c(name = NA_character_, taxid = NA_character_))
  # If a rank occurs more than once, the most specific occurrence is last.
  i <- tail(idx, 1)
  c(name = names[i], taxid = taxids[i])
}

rows <- lapply(seq_len(nrow(lin)), function(i) {
  nm <- strsplit(lin[i, 2], ";", fixed = TRUE)[[1]]
  tx <- strsplit(lin[i, 3], ";", fixed = TRUE)[[1]]
  rk <- strsplit(lin[i, 4], ";", fixed = TRUE)[[1]]
  n <- min(length(nm), length(tx), length(rk))
  nm <- nm[seq_len(n)]; tx <- tx[seq_len(n)]; rk <- rk[seq_len(n)]
  k <- pick_rank(nm, tx, rk, "kingdom")
  p <- pick_rank(nm, tx, rk, "phylum")
  g <- pick_rank(nm, tx, rk, "genus")
  s <- pick_rank(nm, tx, rk, "species")
  data.frame(
    taxid = as.character(lin[i, 1]),
    kindom = ifelse(is.na(k["name"]), NA_character_, paste0(prefix["kingdom"], k["name"])),
    phylum = ifelse(is.na(p["name"]), NA_character_, paste0(prefix["phylum"], p["name"])),
    genus = ifelse(is.na(g["name"]), NA_character_, paste0(prefix["genus"], g["name"])),
    species = ifelse(is.na(s["name"]), NA_character_, paste0(prefix["species"], s["name"])),
    g_taxid = unname(g["taxid"]),
    s_taxid = unname(s["taxid"]),
    stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, rows)
out <- out[!duplicated(out$taxid), , drop = FALSE]
out <- out[order(as.numeric(out$taxid)), , drop = FALSE]

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("[lineage] wrote ", nrow(out), " rows to ", opt$out)
