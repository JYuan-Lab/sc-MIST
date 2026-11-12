#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

option_list <- list(
  make_option("--samples", type = "character", help = "Comma-separated sample IDs"),
  make_option("--hits-dir", dest = "hits_dir", type = "character",
              help = "Directory containing <sample>.cell_line_quantile_hits.tsv"),
  make_option("--kraken-dir", dest = "kraken_dir", type = "character",
              help = "Directory containing first-pass <sample>.kraken.report.txt"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--threshold", type = "double", default = 0.1,
              help = "decontam frequency threshold [default %default]"),
  make_option("--skip-decontam", dest = "skip_decontam", action = "store_true",
              default = FALSE, help = "Write NA p-values instead of running decontam")
)
opt <- parse_args(OptionParser(option_list = option_list))

samples <- trimws(strsplit(opt$samples, ",")[[1]])
samples <- samples[nzchar(samples)]
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

all_hits <- data.frame()
used_samples <- character()
for (s in samples) {
  f <- file.path(opt$hits_dir, paste0(s, ".cell_line_quantile_hits.tsv"))
  if (!file.exists(f)) f <- file.path(opt$hits_dir, s, paste0(s, ".cell_line_quantile_hits.tsv"))
  if (!file.exists(f)) stop("Missing hits file: ", f)
  x <- read.delim(f, header = TRUE, sep = "\t", check.names = FALSE)
  if (nrow(x) == 0) next
  x$run <- s
  all_hits <- rbind(all_hits, x)
  used_samples <- c(used_samples, s)
}

if (nrow(all_hits) == 0) {
  write.table(data.frame(taxid = character(), p = numeric(), contaminant = logical()),
              file.path(opt$outdir, "decontam.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

all_hits$taxid <- as.character(all_hits$taxid)
wide <- all_hits %>%
  select(taxid, run, rpmm) %>%
  pivot_wider(names_from = run, values_from = rpmm, values_fill = 0)

taxa <- as.character(wide$taxid)
count_table <- as.data.frame(t(wide[, -1, drop = FALSE]))
count_table <- count_table[used_samples, , drop = FALSE]
rownames(count_table) <- used_samples
colnames(count_table) <- taxa
count_matrix <- data.matrix(count_table)

concentrations <- sapply(used_samples, function(s) {
  f <- file.path(opt$kraken_dir, paste0(s, ".kraken.report.txt"))
  if (!file.exists(f)) f <- file.path(opt$kraken_dir, s, paste0(s, ".kraken.report.txt"))
  if (!file.exists(f)) stop("Missing first-pass Kraken report: ", f)
  d <- read.delim(f, header = FALSE, sep = "\t", quote = "")
  sum(d$V2[d$V7 %in% c(2, 10239)])
})

if (opt$skip_decontam || !requireNamespace("decontam", quietly = TRUE)) {
  if (!opt$skip_decontam) {
    warning("R package 'decontam' is not installed; writing NA p-values")
  }
  out <- data.frame(taxid = taxa, p = NA_real_, contaminant = FALSE)
} else {
  res <- decontam::isContaminant(
    count_matrix,
    normalize = FALSE,
    conc = concentrations,
    method = "frequency",
    threshold = opt$threshold
  )
  out <- data.frame(taxid = rownames(res), res, check.names = FALSE)
}

write.table(out, file.path(opt$outdir, "decontam.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[09] done: ", nrow(out), " taxa")
