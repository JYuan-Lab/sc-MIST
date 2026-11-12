#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--kraken-report", dest = "kraken_report", type = "character", help = "Kraken2 report"),
  make_option("--mpa-report", dest = "mpa_report", type = "character", help = "MPA-style report"),
  make_option("--kraken-output", dest = "kraken_output", type = "character", help = "Kraken2 per-read output"),
  make_option("--fq1", type = "character", help = "Classified R1 FASTQ"),
  make_option("--fq2", type = "character", default = NULL, help = "Classified R2 FASTQ"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--ntaxid", type = "integer", default = 8000,
              help = "Number of taxids per grep call [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_microbial_taxids <- function(kraken_report, mpa_report) {
  kr <- read.delim(kraken_report, header = FALSE, sep = "\t", quote = "")
  kr$V8 <- trimws(kr$V8)
  kr <- kr[!kr$V8 %in% c("unclassified", "root"), , drop = FALSE]
  kr_key <- gsub("[^[:alnum:]]", "_", kr$V8)

  mpa <- read.delim(mpa_report, header = FALSE, sep = "\t", quote = "")
  micro_lines <- mpa$V1[str_detect(mpa$V1, "(^|\\|)k__(Bacteria|Fungi|Viruses)")]
  if (length(micro_lines) == 0) {
    return(character(0))
  }

  # Match the terminal name of each MPA lineage to the Kraken report instead
  # of relying on row order. This is robust to extra intermediate ranks.
  terminal <- sub(".*\\|", "", micro_lines)
  terminal <- sub("^[a-z]__", "", terminal)
  terminal <- gsub("[^[:alnum:]]", "_", terminal)
  idx <- match(terminal, kr_key)
  as.character(kr$V7[idx[!is.na(idx)]])
}

write_fixed_patterns <- function(values, prefix = "", suffix = "") {
  f <- tempfile()
  writeLines(paste0(prefix, values, suffix), f)
  f
}

extract_microbiome_output <- function(kraken_output, taxids, out_file) {
  if (length(taxids) == 0) {
    file.create(out_file)
    return(invisible(NULL))
  }
  pat <- write_fixed_patterns(taxids, "(taxid ", ")")
  on.exit(unlink(pat), add = TRUE)
  status <- system2("grep", c("-F", "-f", pat, kraken_output), stdout = out_file)
  if (status != 0) {
    stop("grep failed while extracting the Kraken output")
  }
}

extract_microbial_fasta <- function(fq, taxids, sample_name, outdir) {
  out_fa <- file.path(outdir, paste0(sample_name, ".fa"))
  if (length(taxids) == 0) {
    file.create(out_fa)
    return(out_fa)
  }

  pat <- write_fixed_patterns(taxids, "taxid|")
  on.exit(unlink(pat), add = TRUE)

  line_file <- tempfile()
  on.exit(unlink(line_file), add = TRUE)
  status <- system2("grep", c("-F", "-w", "-n", "-f", pat, fq), stdout = line_file)
  if (status != 0) {
    stop("grep failed while finding microbial reads in ", fq)
  }

  h <- sub(":.*$", "", readLines(line_file, warn = FALSE))
  h <- h[nzchar(h)]
  if (length(h) == 0) {
    file.create(out_fa)
    return(out_fa)
  }

  # Keep the FASTQ header and sequence lines only (2 lines per read).
  h <- as.integer(h)
  idx <- sort(unique(c(h, h + 1L)))
  idx_file <- tempfile()
  on.exit(unlink(idx_file), add = TRUE)
  writeLines(as.character(idx), idx_file)

  tmp_fa <- tempfile()
  on.exit(unlink(tmp_fa), add = TRUE)
  status <- system2(
    "awk",
    c(shQuote("NR==FNR{a[$1];next} FNR in a"), shQuote(idx_file), shQuote(fq)),
    stdout = tmp_fa
  )
  if (status != 0) {
    stop("awk failed while extracting reads from ", fq)
  }

  x <- readLines(tmp_fa, warn = FALSE)
  if (length(x) > 0) {
    x[seq(1, length(x), by = 2)] <- sub("^@", ">", x[seq(1, length(x), by = 2)])
  }
  writeLines(x, out_fa)
  out_fa
}

taxids <- unique(read_microbial_taxids(opt$kraken_report, opt$mpa_report))
message("[03] ", opt$sample, ": ", length(taxids), " microbial taxids")

extract_microbiome_output(
  kraken_output = opt$kraken_output,
  taxids = taxids,
  out_file = file.path(opt$outdir, paste0(opt$sample, ".microbiome.output.txt"))
)

invisible(extract_microbial_fasta(opt$fq1, taxids, paste0(opt$sample, "_1"), opt$outdir))
if (!is.null(opt$fq2)) {
  invisible(extract_microbial_fasta(opt$fq2, taxids, paste0(opt$sample, "_2"), opt$outdir))
}

message("[03] done: ", opt$sample)
