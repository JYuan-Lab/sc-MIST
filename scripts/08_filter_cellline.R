#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--barcode-hits", dest = "barcode_hits", type = "character",
              help = "Step 06 <sample>.barcode.kmer.hits.tsv"),
  make_option("--rpmm", type = "character",
              help = "Step 06 <sample>.kraken.report.rpmm.tsv"),
  make_option("--sample-denoising", dest = "sample_denoising", type = "character",
              help = "Step 07 sample_denoising.tsv"),
  make_option("--kraken-report", dest = "kraken_report", type = "character",
              help = "Step 05 Kraken2 report"),
  make_option("--cellline-q90", dest = "cellline_q90", type = "character",
              help = "Cell-line RPMM 90th-percentile table"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--hbv-threshold", dest = "hbv_threshold", type = "double", default = 10.507194,
              help = "HBV (taxid 10407) override threshold [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$barcode_hits, opt$rpmm, opt$sample_denoising,
            opt$kraken_report, opt$cellline_q90)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

q_df <- read.delim(opt$cellline_q90, header = TRUE, sep = "\t", check.names = FALSE)
q_col <- if ("q90" %in% names(q_df)) {
  "q90"
} else if ("threshold_0.9_exclude0" %in% names(q_df)) {
  "threshold_0.9_exclude0"
} else {
  stop("cellline-q90 must contain a 'q90' or 'threshold_0.9_exclude0' column")
}
if (!"taxid" %in% names(q_df)) {
  stop("cellline-q90 must contain a 'taxid' column")
}
q_df <- q_df[, c("taxid", q_col), drop = FALSE]
names(q_df) <- c("taxid", "CLrpmm")
q_df$taxid <- as.character(q_df$taxid)
if (any(q_df$taxid == "10407")) {
  q_df$CLrpmm[q_df$taxid == "10407"] <- opt$hbv_threshold
}

kr <- read.csv(opt$rpmm, header = TRUE, sep = "\t", check.names = FALSE)
kr$taxid <- as.character(kr$taxid)
kr <- subset(kr, rank %in% c("G", "S"))

bh <- read.csv(opt$barcode_hits, header = TRUE, sep = "\t", check.names = FALSE)
bh$taxid <- as.character(bh$taxid)
if ("name" %in% names(bh)) bh$name <- NULL

kr_sub <- kr %>% select(taxid, name, rpmm) %>% distinct(taxid, .keep_all = TRUE)
bh_sub <- bh %>% distinct(taxid, .keep_all = TRUE)
q_sub <- q_df %>% distinct(taxid, .keep_all = TRUE)

merged <- full_join(kr_sub, bh_sub, by = "taxid") %>%
  left_join(q_sub, by = "taxid")
merged$CLrpmm[is.na(merged$CLrpmm)] <- 0

sample_level <- read.csv(opt$sample_denoising, header = TRUE, sep = "\t", check.names = FALSE)
sample_level$taxid <- as.character(sample_level$taxid)
s1 <- left_join(merged, sample_level, by = "taxid")

if ("name.y" %in% names(s1)) {
  s1$name <- ifelse(is.na(s1$name.y), s1$name.x, s1$name.y)
  s1$name.x <- NULL
  s1$name.y <- NULL
} else {
  s1$name <- s1$name.x
  s1$name.x <- NULL
}
if ("rank.y" %in% names(s1)) {
  s1$rank <- ifelse(is.na(s1$rank.y), s1$rank.x, s1$rank.y)
  s1$rank.x <- NULL
  s1$rank.y <- NULL
} else if ("rank.x" %in% names(s1)) {
  s1$rank <- s1$rank.x
  s1$rank.x <- NULL
}

kr_report <- read.delim(opt$kraken_report, header = FALSE, sep = "\t", quote = "")
kr_report$V8 <- trimws(kr_report$V8)
kr_1 <- unique(data.frame(rank = kr_report$V6,
                          taxid = as.character(kr_report$V7),
                          name = kr_report$V8))
s1$rank <- ifelse(s1$taxid %in% kr_1$taxid, kr_1$rank[match(s1$taxid, kr_1$taxid)], s1$rank)
s1$name <- ifelse(s1$taxid %in% kr_1$taxid, kr_1$name[match(s1$taxid, kr_1$taxid)], s1$name)

s1 <- subset(s1, !taxid %in% c("9606", "9605"))
s1$mreads <- ifelse(s1$taxid %in% kr$taxid,
                    kr$reads[match(s1$taxid, kr$taxid)], NA)

write.table(s1, file.path(opt$outdir, paste0(opt$sample, ".cell_line_quantile_hits.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(s1$taxid, file.path(opt$outdir, paste0(opt$sample, ".cell_line_quantile_hits_taxa.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("[08] done: ", opt$sample, " (", nrow(s1), " candidate taxa)")
