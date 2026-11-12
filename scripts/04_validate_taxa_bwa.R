#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

option_list <- list(
  make_option("--sample", type = "character", help = "Sample ID"),
  make_option("--microbiome-output", dest = "microbiome_output", type = "character",
              help = "Step 03 .microbiome.output.txt"),
  make_option("--fa1", type = "character", help = "Step 03 R1 FASTA"),
  make_option("--fa2", type = "character", help = "Step 03 R2 FASTA"),
  make_option("--lineage", type = "character", help = "taxonkit lineage matrix TSV"),
  make_option("--curated-species", dest = "curated_species", type = "character",
              help = "resources/curated_bacteria_species.tsv"),
  make_option("--bwa-library", dest = "bwa_library", type = "character",
              help = "Directory containing {taxid}.fa BWA reference files"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--bwa", type = "character", default = "bwa"),
  make_option("--samtools", type = "character", default = "samtools"),
  make_option("--seqkit", type = "character", default = "seqkit"),
  make_option("--seqtk", type = "character", default = "seqtk"),
  make_option("--bwa-filter", dest = "bwa_filter", type = "character",
              default = file.path(script_dir, "utils", "bwa_filter.awk")),
  make_option("--min-score", dest = "min_score", type = "double", default = 6,
              help = "Minimum curated bacterial score [default %default]"),
  make_option("--min-aln", dest = "min_aln", type = "integer", default = 45,
              help = "Minimum alignment length [default %default]"),
  make_option("--min-identity", dest = "min_identity", type = "double", default = 0.95,
              help = "Minimum sequence identity [default %default]"),
  make_option("--threads", type = "integer", default = 1,
              help = "Number of taxa processed in parallel [default %default]"),
  make_option("--bwa-threads", dest = "bwa_threads", type = "integer", default = 4,
              help = "Threads per BWA process [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$microbiome_output, opt$fa1, opt$fa2, opt$lineage,
            opt$curated_species, opt$bwa_filter)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
for (exe in c(opt$bwa, opt$samtools, opt$seqkit, opt$seqtk)) {
  if (nzchar(Sys.which(exe)) == FALSE) stop("Executable not found: ", exe)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

curated <- fread(opt$curated_species)
if (!all(c("species", "score") %in% names(curated))) {
  stop("curated-species must contain 'species' and 'score' columns")
}
curated_species <- unique(curated[score >= opt$min_score, species])

lineage <- fread(opt$lineage)
required <- c("taxid", "kindom", "phylum", "genus", "species", "g_taxid", "s_taxid")
if (!all(required %in% names(lineage))) {
  stop("lineage matrix is missing columns: ",
       paste(setdiff(required, names(lineage)), collapse = ", "))
}

micro <- fread(opt$microbiome_output, header = FALSE)
if (nrow(micro) == 0) {
  stop("No microbial reads in ", opt$microbiome_output)
}
micro[, taxid := sub(".*\\(taxid ([0-9]+)\\).*", "\\1", V3)]
micro <- micro[grepl("^[0-9]+$", taxid)]

lineage <- lineage[taxid %in% micro$taxid & (genus != "" | species != "")]
lineage[, `:=`(
  taxid = as.character(taxid),
  g_taxid = as.character(g_taxid),
  s_taxid = as.character(s_taxid)
)]

read_map <- merge(
  micro[, .(readid = V2, taxid = as.character(taxid))],
  lineage[, .(taxid, kindom, phylum, genus, species, gtaxid = g_taxid, staxid = s_taxid)],
  by = "taxid"
)

read_map_virus <- read_map[!(kindom == "k__Viruses" & is.na(gtaxid) & is.na(staxid))]
read_map_bact <- read_map[kindom == "k__Bacteria" & (is.na(staxid) | species %in% curated_species)]
read_map <- unique(rbind(read_map_virus, read_map_bact))

setcolorder(read_map, c("taxid", "readid", "kindom", "phylum", "genus", "gtaxid", "staxid", "species"))
map_file <- file.path(opt$outdir, paste0(opt$sample, "_reads_map_lineage.tsv"))
fwrite(read_map, map_file, sep = "\t", quote = FALSE, na = "NA")

if (nrow(read_map) == 0) {
  warning("No reads passed the lineage/curated-species filter")
  fwrite(data.table(taxid = character(), readid = character()),
         file.path(opt$outdir, paste0(opt$sample, "_bwa_pass_reads.tsv")),
         sep = "\t", quote = FALSE)
  file.create(file.path(opt$outdir, paste0(opt$sample, "_species_1.fq")))
  file.create(file.path(opt$outdir, paste0(opt$sample, "_species_2.fq")))
  fwrite(data.table(taxid = character()),
         file.path(opt$outdir, paste0(opt$sample, "_bwa_black_species_id.tsv")),
         sep = "\t", quote = FALSE, col.names = FALSE)
  quit(save = "no", status = 0)
}

taxa <- unique(c(read_map$gtaxid[!is.na(read_map$gtaxid)],
                 read_map$staxid[!is.na(read_map$staxid)]))

validate_one <- function(ta) {
  ref <- file.path(opt$bwa_library, paste0(ta, ".fa"))
  if (!file.exists(ref)) {
    warning("Missing BWA reference: ", ref)
    return(data.table(taxid = ta, readid = character()))
  }

  ids <- if (ta %in% read_map$gtaxid) {
    read_map[gtaxid == ta, readid]
  } else {
    read_map[staxid == ta, readid]
  }
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) == 0) {
    return(data.table(taxid = ta, readid = character()))
  }

  id_file <- tempfile()
  sub_fq <- tempfile(fileext = ".fq")
  on.exit(unlink(c(id_file, sub_fq)), add = TRUE)
  writeLines(ids, id_file)

  cmd_extract <- paste(
    shQuote(opt$seqkit), "grep -j 1 -f", shQuote(id_file), shQuote(opt$fa2),
    "|", shQuote(opt$seqkit), "replace -p", shQuote(" .*"), "-w 0",
    "|", shQuote(opt$seqtk), "seq -F", shQuote("I"), "-",
    ">", shQuote(sub_fq)
  )
  status <- system(cmd_extract)
  if (status != 0 || !file.exists(sub_fq) || file.info(sub_fq)$size == 0) {
    return(data.table(taxid = ta, readid = character()))
  }

  cmd_align <- paste(
    shQuote(opt$bwa), "mem -t", opt$bwa_threads, shQuote(ref), shQuote(sub_fq), "2>/dev/null",
    "|", shQuote(opt$samtools), "view -F 0x900 -",
    "| awk ", paste0("-v min_len=", opt$min_aln), " ",
    paste0("-v min_id=", opt$min_identity),
    " -f", shQuote(opt$bwa_filter)
  )
  pass <- suppressWarnings(system(cmd_align, intern = TRUE))
  pass <- unique(pass[nzchar(pass)])
  data.table(taxid = ta, readid = pass)
}

message("[04] ", opt$sample, ": validating ", length(taxa), " taxa")
res <- mclapply(taxa, validate_one, mc.cores = max(1L, opt$threads))
pass_dt <- unique(rbindlist(res, fill = TRUE))
pass_dt <- pass_dt[nzchar(readid)]

fwrite(pass_dt,
       file.path(opt$outdir, paste0(opt$sample, "_bwa_pass_reads.tsv")),
       sep = "\t", quote = FALSE)

pass_ids <- unique(pass_dt$readid)
extract_pass <- function(fa, out_fq) {
  if (length(pass_ids) == 0) {
    file.create(out_fq)
    return(invisible(NULL))
  }
  id_file <- tempfile()
  on.exit(unlink(id_file), add = TRUE)
  writeLines(pass_ids, id_file)
  cmd <- paste(
    shQuote(opt$seqkit), "grep -j 1 -f", shQuote(id_file), shQuote(fa),
    "|", shQuote(opt$seqkit), "replace -p", shQuote(" .*"), "-w 0",
    "|", shQuote(opt$seqtk), "seq -F", shQuote("I"), "-",
    ">", shQuote(out_fq)
  )
  status <- system(cmd)
  if (status != 0) stop("Failed to extract BWA-passing reads from ", fa)
  invisible(NULL)
}

extract_pass(opt$fa1, file.path(opt$outdir, paste0(opt$sample, "_species_1.fq")))
extract_pass(opt$fa2, file.path(opt$outdir, paste0(opt$sample, "_species_2.fq")))

species_taxids <- unique(read_map$staxid[!is.na(read_map$staxid)])
validated_taxids <- unique(pass_dt$taxid)
black_species <- setdiff(species_taxids, validated_taxids)
fwrite(data.table(taxid = black_species),
       file.path(opt$outdir, paste0(opt$sample, "_bwa_black_species_id.tsv")),
       sep = "\t", quote = FALSE, col.names = FALSE)

message("[04] done: ", opt$sample,
        " (", length(unique(pass_dt$taxid)), " validated taxa; ",
        length(black_species), " black species)")
