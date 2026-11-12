#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--fa1", type = "character", help = "Step 05 microbial R1 FASTA"),
  make_option("--fa2", type = "character", help = "Step 05 microbial R2 FASTA"),
  make_option("--taxa", type = "character", help = "Step 08 taxa list"),
  make_option("--kraken-report", dest = "kraken_report", type = "character",
              help = "Step 05 Kraken2 report"),
  make_option("--mpa-report", dest = "mpa_report", type = "character",
              help = "Step 05 MPA report"),
  make_option("--hits", type = "character", help = "Step 08 cell_line_quantile_hits.tsv"),
  make_option("--decontam", type = "character", help = "Step 09 decontam.tsv"),
  make_option("--bwa-black", dest = "bwa_black", type = "character", default = NULL,
              help = "Step 04 <sample>_bwa_black_species_id.tsv"),
  make_option("--black-genus", dest = "black_genus", type = "character",
              help = "resources/black_list.txt"),
  make_option("--mgnify", type = "character", help = "resources/mgnify_black.txt"),
  make_option("--ccle-blacklist", dest = "ccle_blacklist", type = "character",
              help = "resources/ccle_blacklist_taxids.txt"),
  make_option("--cellranger-outs", dest = "cellranger_outs", type = "character", default = NULL,
              help = "Cell Ranger outs directory (optional, for barcode filtering)"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--cb-len", dest = "cb_len", type = "integer", default = 16),
  make_option("--umi-len", dest = "umi_len", type = "integer", default = 10),
  make_option("--min-corr", dest = "min_corr", type = "double", default = 0,
              help = "Minimum r/r1/r2/r3 [default %default]"),
  make_option("--decontam-threshold", dest = "decontam_threshold", type = "double", default = 0.1,
              help = "decontam p-value threshold [default %default]"),
  make_option("--mgnify-min-studies", dest = "mgnify_min_studies", type = "integer", default = 4,
              help = "Flag MGnify taxa seen in > this many studies [default %default]"),
  make_option("--no-cellranger-filter", dest = "no_cellranger_filter", action = "store_true",
              default = FALSE, help = "Do not restrict output to Cell Ranger barcodes"),
  make_option("--nfilter", type = "integer", default = 130,
              help = "Remove reads with >= this many identical bases [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$fa1, opt$fa2, opt$taxa, opt$kraken_report,
            opt$mpa_report, opt$hits, opt$decontam,
            opt$black_genus, opt$mgnify, opt$ccle_blacklist)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (file.info(opt$fa1)$size == 0 || file.info(opt$fa2)$size == 0) {
  empty_counts <- data.table(barcode = character(), counts = integer(), taxid = character(),
                             rank = character(), kingdom = character(), phylum = character(),
                             class = character(), order = character(), family = character(),
                             genus = character(), species = character())
  fwrite(empty_counts, file.path(opt$outdir, paste0(opt$sample, ".counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(empty_counts, file.path(opt$outdir, paste0(opt$sample, ".filtered_counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(data.table(taxid = character(), rank = character(), name = character(),
                    genus = character(), species = character(), umi = integer()),
         file.path(opt$outdir, paste0(opt$sample, ".taxon_abundance.tsv")),
         sep = "\t", quote = FALSE)
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

ok <- nchar(seq1) >= (opt$cb_len + opt$umi_len)
seq1 <- seq1[ok]; h1 <- h1[ok]
barcode <- substr(seq1, 1, opt$cb_len)
umi <- substr(seq1, opt$cb_len + 1, opt$cb_len + opt$umi_len)
taxid <- sub(".*taxid\\|", "", h1)
id <- sub("\\s.*$", "", h1)

keep_pair <- id %in% sub("\\s.*$", "", h2)
id <- id[keep_pair]; barcode <- barcode[keep_pair]; umi <- umi[keep_pair]
taxid <- taxid[keep_pair]; seq1 <- seq1[keep_pair]
seq2 <- seq2[match(id, sub("\\s.*$", "", h2))]

if (opt$nfilter > 0 && length(seq1) > 0) {
  keep <- !(
    grepl(paste0("A{", opt$nfilter, "}"), seq1) |
    grepl(paste0("C{", opt$nfilter, "}"), seq1) |
    grepl(paste0("G{", opt$nfilter, "}"), seq1) |
    grepl(paste0("T{", opt$nfilter, "}"), seq1)
  )
  id <- id[keep]; barcode <- barcode[keep]; umi <- umi[keep]
  taxid <- taxid[keep]; seq1 <- seq1[keep]; seq2 <- seq2[keep]
}

taxa <- if (file.info(opt$taxa)$size > 0) {
  fread(opt$taxa, header = FALSE)$V1 %>% as.character()
} else {
  character()
}

write_empty <- function() {
  empty_counts <- data.table(barcode = character(), counts = integer(), taxid = character(),
                             rank = character(), kingdom = character(), phylum = character(),
                             class = character(), order = character(), family = character(),
                             genus = character(), species = character())
  fwrite(empty_counts, file.path(opt$outdir, paste0(opt$sample, ".counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(empty_counts, file.path(opt$outdir, paste0(opt$sample, ".filtered_counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(data.table(taxid = character(), rank = character(), name = character(),
                    genus = character(), species = character(), umi = integer()),
         file.path(opt$outdir, paste0(opt$sample, ".taxon_abundance.tsv")),
         sep = "\t", quote = FALSE)
}

if (length(taxa) == 0) {
  write_empty()
  quit(save = "no", status = 0)
}

kr <- fread(opt$kraken_report, header = FALSE)
kr[, V8 := trimws(V8)]
kr <- kr[!V8 %in% c("unclassified", "root")]
kr[, V8 := str_replace_all(V8, "[^[:alnum:]]", "_")]
kr[, V7 := as.character(V7)]
kr.names <- str_squish(str_replace_all(str_remove(kr$V8, ".*__"), "_", " "))

mpa <- fread(opt$mpa_report, header = FALSE)
mpa[, taxid := NA_character_]
if (nrow(mpa) > 1) {
  for (i in seq.int(2, nrow(mpa))) {
    t_names <- mpa$V1[i] %>% as.character() %>% strsplit("\\|") %>% unlist() %>%
      str_remove(".*__") %>% str_replace_all("[^[:alnum:]]", "_")
    mpa$taxid[i] <- paste0("*", paste(kr$V7[match(t_names, kr$V8)], collapse = "*"), "*")
  }
}

x <- paste0("\\*", taxa, "\\*", collapse = "|")
ii <- str_which(mpa$taxid, x)
x2 <- paste0("\\*", taxa, "\\*.*")
ii <- sapply(x2, function(y) str_extract(mpa$taxid[ii], y) %>% strsplit("\\*") %>%
               unlist() %>% unique()) %>% unlist() %>% unique()

i <- which(taxid %in% ii)
barcode <- barcode[i]; umi <- umi[i]; taxid <- taxid[i]

bc <- data.frame(barcode = barcode, umi = umi, taxid = taxid, stringsAsFactors = FALSE)
s.bc <- unique(bc)
s.bc$umi <- 1L
s.bc <- s.bc %>% group_by(barcode, taxid) %>% summarize(umi = sum(umi), .groups = "keep")

dummy <- data.table(barcode = character(), counts = integer(), taxid = character())
if (nrow(s.bc) == 0) {
  fwrite(dummy, file.path(opt$outdir, paste0(opt$sample, ".counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(dummy, file.path(opt$outdir, paste0(opt$sample, ".filtered_counts.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(data.table(taxid = character(), rank = character(), name = character(),
                    genus = character(), species = character(), umi = integer()),
         file.path(opt$outdir, paste0(opt$sample, ".taxon_abundance.tsv")),
         sep = "\t", quote = FALSE)
  quit(save = "no", status = 0)
}

s.bc$barcode_factor <- factor(s.bc$barcode)
s.bc$taxid_factor <- factor(s.bc$taxid)
s.mat <- sparseMatrix(i = as.integer(s.bc$barcode_factor),
                      j = as.integer(s.bc$taxid_factor),
                      x = s.bc$umi)
colnames(s.mat) <- levels(s.bc$taxid_factor)
rownames(s.mat) <- levels(s.bc$barcode_factor)
s.mat <- t(s.mat)

count_rows <- list()
counter <- 0L
for (i in seq_len(nrow(s.mat))) {
  row_taxid <- rownames(s.mat)[i]
  hit <- which(kr$V7 == row_taxid)
  if (length(hit) == 0) next
  tax <- mpa[hit[1], 1] %>% as.character() %>% strsplit("\\|") %>% unlist() %>%
    data.frame() %>% separate(1, c("rank", "name"), sep = "__") %>%
    subset(rank %in% c("k", "p", "c", "o", "f", "g", "s")) %>%
    mutate(name = str_replace_all(name, "[^[:alnum:]]", " ") %>% str_squish())
  if (nrow(tax) == 0) next
  tax$id <- vapply(tax$name, function(nm) kr$V7[kr.names == nm][1], character(1))
  row <- s.mat[i, ]
  row <- row[row > 0]
  for (j in seq_len(nrow(tax))) {
    tax2 <- tax[1:j, ]
    counter <- counter + 1L
    count_rows[[counter]] <- data.frame(
      barcode = names(row), counts = as.integer(row), taxid = tax2$id[j],
      rank = tax2$rank[j],
      kingdom = tax2$name[tax2$rank == "k"][1],
      phylum = tax2$name[tax2$rank == "p"][1],
      class = tax2$name[tax2$rank == "c"][1],
      order = tax2$name[tax2$rank == "o"][1],
      family = tax2$name[tax2$rank == "f"][1],
      genus = tax2$name[tax2$rank == "g"][1],
      species = tax2$name[tax2$rank == "s"][1]
    )
  }
}

counts <- rbindlist(count_rows, fill = TRUE)
if (nrow(counts) == 0) {
  write_empty()
  quit(save = "no", status = 0)
}
counts <- counts %>%
  group_by(barcode, taxid, rank, kingdom, phylum, class, order, family, genus, species) %>%
  summarize(counts = sum(counts), .groups = "drop") %>%
  arrange(desc(counts))
counts <- as.data.table(counts)
fwrite(counts, file.path(opt$outdir, paste0(opt$sample, ".counts.tsv")),
       sep = "\t", quote = FALSE)

lineage_tab <- unique(counts[, c("taxid", "genus", "species")])

hits <- fread(opt$hits)
hits[, taxid := as.character(taxid)]
hits[, rank := toupper(rank)]
metrics <- merge(hits, lineage_tab, by = "taxid", all.x = TRUE)

dc <- fread(opt$decontam)
dc[, taxid := as.character(taxid)]
metrics <- merge(metrics, dc[, .(taxid, decontam_p = p)], by = "taxid", all.x = TRUE)

black1 <- read.delim(opt$black_genus, header = TRUE, sep = "\t", check.names = FALSE)
mgnify <- read.csv(opt$mgnify, header = TRUE, sep = ",", check.names = FALSE)
mgnify <- subset(mgnify, count_n > opt$mgnify_min_studies)
black_names <- c(trimws(black1$genus), trimws(mgnify$taxonomy))
metrics$black <- ifelse(metrics$genus %in% black_names | metrics$species %in% black_names,
                        "YES", "NO")

bwa_black <- character()
if (!is.null(opt$bwa_black) && file.exists(opt$bwa_black) && file.info(opt$bwa_black)$size > 0) {
  bwa_black <- fread(opt$bwa_black, header = FALSE)$V1 %>% as.character()
}
ccle_black <- fread(opt$ccle_blacklist, header = FALSE)$V1 %>% as.character()

metrics$pass <- with(metrics,
  rank %in% c("G", "S") &
  !is.na(r) & r > opt$min_corr &
  !is.na(r1) & r1 > opt$min_corr &
  !is.na(r2) & r2 > opt$min_corr &
  !is.na(r3) & r3 > opt$min_corr &
  !is.na(rpmm) & !is.na(CLrpmm) & rpmm > CLrpmm &
  black == "NO" &
  (is.na(decontam_p) | decontam_p > opt$decontam_threshold) &
  !taxid %in% bwa_black & !taxid %in% ccle_black
)

fwrite(metrics, file.path(opt$outdir, paste0(opt$sample, ".taxon_metrics.tsv")),
       sep = "\t", quote = FALSE)

keep_taxa <- metrics$taxid[metrics$pass]
filtered <- counts[taxid %in% keep_taxa & rank %in% c("g", "s")]
filtered <- merge(filtered, metrics[, .(taxid, name, r, r1, r2, r3, rpmm, CLrpmm, decontam_p)],
                  by = "taxid", all.x = TRUE)

if (!is.null(opt$cellranger_outs) && dir.exists(opt$cellranger_outs) && !opt$no_cellranger_filter) {
  bc_file <- file.path(opt$cellranger_outs, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  if (!file.exists(bc_file)) bc_file <- file.path(opt$cellranger_outs, "barcodes.tsv.gz")
  if (file.exists(bc_file)) {
    barcodes <- readLines(gzfile(bc_file))
    bc1 <- paste0(filtered$barcode, "-1")
    filtered$barcode_cellranger <- ifelse(bc1 %in% barcodes, bc1, filtered$barcode)
    filtered$in_cellranger <- filtered$barcode_cellranger %in% barcodes
    filtered <- filtered[filtered$in_cellranger, ]
  } else {
    warning("No barcodes.tsv.gz found under ", opt$cellranger_outs)
  }
}

fwrite(filtered, file.path(opt$outdir, paste0(opt$sample, ".filtered_counts.tsv")),
       sep = "\t", quote = FALSE)

abundance <- if (nrow(filtered) == 0) {
  data.table(taxid = character(), rank = character(), name = character(),
             genus = character(), species = character(), umi = integer())
} else {
  filtered[, .(umi = sum(counts)), by = .(taxid, rank, name, genus, species,
                                          r, r1, r2, r3, rpmm, CLrpmm, decontam_p)]
}
fwrite(abundance, file.path(opt$outdir, paste0(opt$sample, ".taxon_abundance.tsv")),
       sep = "\t", quote = FALSE)

qc <- data.frame(
  sample = opt$sample,
  microbial_umi = sum(abundance$umi),
  taxa_detected = nrow(abundance),
  barcodes_with_microbes = if (nrow(filtered) > 0) length(unique(filtered$barcode)) else 0
)
write.table(qc, file.path(opt$outdir, paste0(opt$sample, ".qc.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("[10] done: ", opt$sample, " (", nrow(abundance), " retained taxa; ",
        sum(abundance$umi), " UMI)")
