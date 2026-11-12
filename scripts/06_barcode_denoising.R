#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--fa1", type = "character", help = "Step 05 microbial R1 FASTA"),
  make_option("--fa2", type = "character", help = "Step 05 microbial R2 FASTA"),
  make_option("--kraken-report", dest = "kraken_report", type = "character",
              help = "Step 05 Kraken2 report"),
  make_option("--mpa-report", dest = "mpa_report", type = "character",
              help = "Step 05 MPA report"),
  make_option("--microbiome-output", dest = "microbiome_output", type = "character",
              help = "Step 05 microbiome output"),
  make_option("--first-kraken-report", dest = "first_kraken_report", type = "character",
              help = "Step 02 Kraken2 report (for RPMM denominator)"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--cb-len", dest = "cb_len", type = "integer", default = 16),
  make_option("--umi-len", dest = "umi_len", type = "integer", default = 10),
  make_option("--min-frac", dest = "min_frac", type = "double", default = 0.5,
              help = "Minimum fraction of k-mers assigned to the taxon [default %default]"),
  make_option("--kmer-len", dest = "kmer_len", type = "integer", default = 35),
  make_option("--max-barcodes", dest = "max_barcodes", type = "integer", default = 1000),
  make_option("--host", type = "character", default = "9606")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$fa1, opt$fa2, opt$kraken_report, opt$mpa_report,
            opt$microbiome_output, opt$first_kraken_report)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

sckmer_out <- file.path(opt$outdir, paste0(opt$sample, ".sckmer.txt"))
barcode_out <- file.path(opt$outdir, paste0(opt$sample, ".barcode.kmer.hits.tsv"))
rpmm_out <- file.path(opt$outdir, paste0(opt$sample, ".kraken.report.rpmm.tsv"))

if (file.info(opt$fa1)$size == 0 || file.info(opt$fa2)$size == 0) {
  fwrite(data.table(barcode = character(), taxid = character(), kmer = integer(), uniq = integer()),
         sckmer_out, sep = "\t", quote = FALSE)
  fwrite(data.table(taxid = character(), r = numeric(), p = numeric(), name = character()),
         barcode_out, sep = "\t", quote = FALSE)
  fwrite(data.table(study = character(), sample = character(), rank = character(),
                    taxid = character(), name = character(), reads = integer(),
                    min = integer(), uniq = integer(), rpm = numeric(), rpmm = numeric(),
                    mreads = integer()),
         rpmm_out, sep = "\t", quote = FALSE)
  message("[06] ", opt$sample, ": no reads after host filtering; wrote empty tables")
  quit(save = "no", status = 0)
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else "."
}
source(file.path(script_dir(), "utils", "read_fasta.R"))

reads1 <- read_fasta(opt$fa1)
reads2 <- read_fasta(opt$fa2)
seq1 <- as.character(reads1)
seq2 <- as.character(reads2)
h1 <- names(reads1)
h2 <- names(reads2)

barcode <- substr(seq1, 1, opt$cb_len)
umi <- substr(seq1, opt$cb_len + 1, opt$cb_len + opt$umi_len)
taxid <- sub(".*taxid\\|", "", h1)
id <- sub("\\s.*$", "", h1)
id2 <- sub("\\s.*$", "", h2)

ok <- nchar(seq1) >= (opt$cb_len + opt$umi_len) & grepl("^[0-9]+$", taxid)
barcode <- barcode[ok]
umi <- umi[ok]
taxid <- taxid[ok]
id <- id[ok]
seq1 <- seq1[ok]

keep <- id %in% id2
id <- id[keep]
barcode <- barcode[keep]
umi <- umi[keep]
taxid <- taxid[keep]
seq1 <- seq1[keep]
idx2 <- match(id, id2)
seq2 <- seq2[idx2]

message("[06] ", opt$sample, ": ", length(id), " paired microbial reads")

kr <- fread(opt$kraken_report, header = FALSE)
kr[, V8 := trimws(V8)]
kr <- kr[!V8 %in% c("unclassified", "root")]
kr[, V8 := str_replace_all(V8, "[^[:alnum:]]", "_")]
kr[, V7 := as.character(V7)]

mpa <- fread(opt$mpa_report, header = FALSE)
mpa[, taxid := NA_character_]
if (nrow(mpa) > 1) {
  for (i in seq.int(2, nrow(mpa))) {
    t_names <- mpa$V1[i] %>%
      as.character() %>%
      strsplit("\\|") %>%
      unlist() %>%
      str_remove(".*__") %>%
      str_replace_all("[^[:alnum:]]", "_")
    mpa$taxid[i] <- paste0("*", paste(kr$V7[match(t_names, kr$V8)], collapse = "*"), "*")
  }
}

micro <- fread(opt$microbiome_output, header = FALSE)
micro <- micro %>%
  select(-V1) %>%
  separate(V3, into = c("name", "taxid"), sep = "\\(taxid") %>%
  mutate(taxid = str_remove(taxid, "\\)") %>% trimws(),
         name = trimws(name))

tx <- kr$V7[kr$V6 %in% c("S", "G")] %>% setdiff(opt$host)
tx <- micro$taxid[micro$taxid %in% tx] %>% unique()

barcode_kmer <- list()
counter <- 0L

for (ta in tx) {
  counter <- counter + 1L
  cat("\r", "taxa processed:", round(counter / length(tx) * 100, 2), "%")

  lin <- str_subset(mpa$taxid, paste0("\\*", ta, "\\*")) %>%
    str_extract(paste0("\\*", ta, "\\*.*")) %>%
    str_remove("^\\*") %>%
    str_remove("\\*$") %>%
    str_split("\\*") %>%
    unlist() %>%
    as.numeric() %>%
    unique()

  full.lin <- str_subset(mpa$taxid, paste0("\\*", ta, "\\*")) %>%
    str_remove("^\\*") %>%
    str_remove("\\*$") %>%
    str_split("\\*") %>%
    unlist() %>%
    as.numeric() %>%
    unique()

  out <- subset(micro, taxid %in% as.character(lin)) %>%
    separate(V5, into = c("r1", "r2"), sep = "\\|\\:\\|")
  out$r1[str_which(out$r1, paste0(" ", opt$host, ":"))] <- NA
  out$r2[str_which(out$r2, paste0(" ", opt$host, ":"))] <- NA
  out$r1[out$r1 == ""] <- NA
  out$r2[out$r2 == ""] <- NA
  out <- subset(out, !is.na(r1) | !is.na(r2))
  out$r1 <- trimws(out$r1)
  out$r2 <- trimws(out$r2)

  if (nrow(out) == 0) next
  i <- which(id %in% out$V2)
  if (length(i) == 0) next

  tax.df <- list()
  n <- if (nrow(out) > opt$max_barcodes) sample(nrow(out), opt$max_barcodes) else seq_len(nrow(out))
  counter2 <- 0L

  for (j in n) {
    for (mate in c("r1", "r2")) {
      r <- data.frame(pos = out[[mate]][j] %>% strsplit("\\s") %>% unlist()) %>%
        separate(pos, into = c("taxid", "nkmer"), sep = ":", convert = TRUE)
      r <- r %>% filter(!is.na(taxid), !is.na(nkmer))
      if (nrow(r) == 0) next

      r <- r %>%
        mutate(fkmer = nkmer / sum(nkmer),
               nt_start = cumsum(nkmer) - nkmer + 1,
               nt_end = cumsum(nkmer) + opt$kmer_len - 1) %>%
        mutate(nt_len = nt_end - nt_start + 1) %>%
        subset(taxid %in% full.lin)

      if (sum(r$fkmer, na.rm = TRUE) < opt$min_frac) next
      if (nrow(r) == 0) next
      counter2 <- counter2 + 1L

      kmer <- character()
      seq_mate <- unname(if (mate == "r1") seq1[i[j]] else seq2[i[j]])
      for (k in seq_len(nrow(r))) {
        if (is.na(r$nkmer[k]) || r$nkmer[k] < 1) next
        for (m in seq_len(r$nkmer[k])) {
          start <- r$nt_start[k] + m - 1
          end <- start + opt$kmer_len - 2
          if (!is.na(end) && end <= nchar(seq_mate)) {
            kmer <- c(kmer, unname(substr(seq_mate, start, end)))
          }
        }
      }
      if (length(kmer) > 0) {
        tax.df[[counter2]] <- data.frame(
          barcode = unname(barcode[i[j]]), taxid = ta, k = unname(kmer),
          n = as.numeric(sum(r$nt_len[r$taxid %in% lin], na.rm = TRUE))
        )
      }
    }
  }

  if (length(tax.df) == 0) next
  barcode_kmer[[counter]] <- bind_rows(tax.df) %>%
    subset(!is.na(k)) %>%
    group_by(barcode, taxid) %>%
    summarize(kmer = n(), uniq = n_distinct(k), .groups = "keep")
}
cat("\n")

if (length(barcode_kmer) == 0) {
  bk <- data.table(barcode = character(), taxid = character(), kmer = integer(), uniq = integer())
} else {
  bk <- rbindlist(barcode_kmer, fill = TRUE)
}
fwrite(bk, sckmer_out, sep = "\t", quote = FALSE)

cor_tab <- bk %>%
  filter(kmer > 1) %>%
  group_by(taxid) %>%
  summarize(
    nn = n(),
    r = if (n() > 3 && sum(is.finite(kmer) & is.finite(uniq)) >= 2)
      suppressWarnings(as.numeric(cor.test(kmer, uniq, method = "spearman")$estimate)) else NA_real_,
    p = if (n() > 3 && sum(is.finite(kmer) & is.finite(uniq)) >= 2)
      suppressWarnings(cor.test(kmer, uniq, method = "spearman")$p.value) else NA_real_,
    .groups = "keep"
  ) %>%
  mutate(p = p.adjust(p)) %>%
  select(-nn)
cor_tab$name <- kr$V8[match(cor_tab$taxid, kr$V7)]
fwrite(cor_tab, barcode_out, sep = "\t", quote = FALSE)

first <- fread(opt$first_kraken_report, header = FALSE)
first[, V8 := trimws(V8)]
n_microbiome_reads <- sum(first$V2[first$V8 %in% c("Bacteria", "Fungi", "Viruses")])
if (n_microbiome_reads == 0) n_microbiome_reads <- NA_real_

total_reads <- first$V2[first$V8 == "root"]
if (length(total_reads) == 0) total_reads <- sum(first$V2[1:min(2, nrow(first))])

rpmm <- kr[, .(
  study = NA_character_, sample = opt$sample, rank = V6,
  taxid = V7, name = V8, reads = V2, min = V4, uniq = V5,
  rpm = V2 / total_reads * 1e6,
  rpmm = V2 / n_microbiome_reads * 1e6,
  mreads = V2
)]
rpmm <- rpmm[reads >= 2 & uniq >= 2]
fwrite(rpmm, rpmm_out, sep = "\t", quote = FALSE)

message("[06] done: ", opt$sample)
