#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

option_list <- list(
  make_option("--samples", type = "character",
              help = "Comma-separated sample IDs (at least 3)"),
  make_option("--rpmm-dir", dest = "rpmm_dir", type = "character",
              help = "Directory containing <sample>.kraken.report.rpmm.tsv"),
  make_option("--outdir", type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

samples <- trimws(strsplit(opt$samples, ",")[[1]])
samples <- samples[nzchar(samples)]
if (length(samples) < 3) {
  warning("Sample-level correlations are under-powered with fewer than 3 samples")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

all <- data.frame()
for (s in samples) {
  f <- file.path(opt$rpmm_dir, paste0(s, ".kraken.report.rpmm.tsv"))
  if (!file.exists(f)) f <- file.path(opt$rpmm_dir, s, paste0(s, ".kraken.report.rpmm.tsv"))
  if (!file.exists(f)) stop("Missing RPMM file: ", f)
  x <- read.csv(f, header = TRUE, sep = "\t", check.names = FALSE)
  if (nrow(x) == 0) next
  x$study <- ""
  x$sample <- s
  all <- rbind(all, x)
}

if (nrow(all) == 0) {
  out <- data.frame(taxid = character(), r1 = numeric(), r2 = numeric(), r3 = numeric(),
                    p1 = numeric(), p2 = numeric(), p3 = numeric(),
                    rank = character(), name = character(), nn = integer())
  write.table(out, file.path(opt$outdir, "sample_denoising.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

all$taxid <- as.character(all$taxid)
all$rank <- as.character(all$rank)

safe_cor <- function(x, y) {
  tryCatch(suppressWarnings(cor(x, y, method = "spearman")), error = function(e) NA_real_)
}
safe_p <- function(x, y) {
  tryCatch(suppressWarnings(cor.test(x, y, method = "spearman")$p.value), error = function(e) NA_real_)
}

kr3 <- all %>%
  group_by(taxid) %>%
  filter(n() > 2, rank %in% c("G", "S")) %>%
  summarize(
    r1 = safe_cor(min, uniq),
    r2 = safe_cor(min, reads),
    r3 = safe_cor(reads, uniq),
    p1 = safe_p(min, uniq),
    p2 = safe_p(min, reads),
    p3 = safe_p(reads, uniq),
    .groups = "keep"
  )

kr_lt3 <- all %>%
  group_by(taxid) %>%
  filter(n() < 3, rank %in% c("G", "S")) %>%
  summarize(r1 = NA_real_, r2 = NA_real_, r3 = NA_real_,
            p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, .groups = "keep")

c2 <- bind_rows(kr3, kr_lt3)

meta <- all %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  select(rank, taxid, name, nn) %>%
  distinct()

out <- c2 %>%
  left_join(meta, by = "taxid") %>%
  filter(rank %in% c("G", "S"))

write.table(out, file.path(opt$outdir, "sample_denoising.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[07] done: ", nrow(out), " taxa across ", length(samples), " samples")
