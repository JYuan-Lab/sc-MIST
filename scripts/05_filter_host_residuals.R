#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--fq1", type = "character", help = "R1 FASTQ after BWA validation"),
  make_option("--fq2", type = "character", help = "R2 FASTQ after BWA validation"),
  make_option("--reference", type = "character",
              help = "Combined HLA + human rRNA FASTA (BWA-indexed)"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--bwa", type = "character", default = "bwa"),
  make_option("--samtools", type = "character", default = "samtools"),
  make_option("--seqkit", type = "character", default = "seqkit"),
  make_option("--threads", type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$fq1, opt$fq2, opt$reference)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
if (!file.exists(paste0(opt$reference, ".bwt"))) {
  stop("BWA index not found for ", opt$reference,
       ". Run `bwa index ", opt$reference, "` first.")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

keep_ids <- tempfile()
sam_ids <- tempfile()
on.exit(unlink(c(keep_ids, sam_ids)), add = TRUE)

cmd_align <- paste(
  shQuote(opt$bwa), "mem -t", opt$threads, shQuote(opt$reference), shQuote(opt$fq2),
  "2>/dev/null |", shQuote(opt$samtools), "view -f 4 - | cut -f1 | sort -u",
  ">", shQuote(keep_ids)
)
status <- system(cmd_align)
if (status != 0) stop("BWA host-residual filtering failed")

clean1 <- file.path(opt$outdir, paste0(opt$sample, "_clean_1.fq"))
clean2 <- file.path(opt$outdir, paste0(opt$sample, "_clean_2.fq"))

if (file.info(keep_ids)$size == 0) {
  file.create(clean1, clean2)
} else {
  cmd1 <- paste(shQuote(opt$seqkit), "grep -j", opt$threads, "-f", shQuote(keep_ids),
                shQuote(opt$fq1), ">", shQuote(clean1))
  cmd2 <- paste(shQuote(opt$seqkit), "grep -j", opt$threads, "-f", shQuote(keep_ids),
                shQuote(opt$fq2), ">", shQuote(clean2))
  if (system(cmd1) != 0) stop("Failed to extract clean R1")
  if (system(cmd2) != 0) stop("Failed to extract clean R2")
}

n_before <- as.integer(system(paste("grep -c '^@'", shQuote(opt$fq1)), intern = TRUE))
n_after <- as.integer(system(paste("grep -c '^@'", shQuote(clean1)), intern = TRUE))
stats <- data.frame(sample = opt$sample, reads_before = n_before,
                    reads_after = n_after, reads_removed = n_before - n_after)
write.table(stats, file.path(opt$outdir, paste0(opt$sample, "_host_filter_stats.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("[05] done: ", opt$sample, " (", n_before, " -> ", n_after, " reads)")
