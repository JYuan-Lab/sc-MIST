#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--samples", type = "character", help = "Comma-separated sample IDs"),
  make_option("--abundance-dir", dest = "abundance_dir", type = "character",
              help = "Directory containing <sample>.taxon_abundance.tsv"),
  make_option("--counts-dir", dest = "counts_dir", type = "character",
              help = "Directory containing <sample>.filtered_counts.tsv"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--lineage", type = "character", default = NULL,
              help = "Seven-column lineage matrix used to annotate kingdoms"),
  make_option("--ncbi-virus-host", dest = "ncbi_virus_host", type = "character",
              default = NULL, help = "NCBI Virus nucleotide host metadata TSV"),
  make_option("--ncbi-virus-assembly-host", dest = "ncbi_virus_assembly_host",
              type = "character", default = NULL,
              help = "Optional NCBI Virus assembly host metadata table"),
  make_option("--hvd", type = "character", default = NULL,
              help = "humanvirusbase_Virus_Organ_pairs.csv"),
  make_option("--hvidb", type = "character", default = NULL,
              help = "HVIDB_PPIs.txt"),
  make_option("--virus-min-score", dest = "virus_min_score", type = "integer",
              default = 2, help = "Keep viruses supported by >= this many sources [default %default]"),
  make_option("--keep-other-kingdoms", dest = "keep_other_kingdoms",
              action = "store_true", default = FALSE,
              help = "Also keep non-bacterial, non-viral taxa (e.g. fungi)"),
  make_option("--no-virus-filter", dest = "no_virus_filter",
              action = "store_true", default = FALSE,
              help = "Write matrices without the viral host-evidence filter")
)
opt <- parse_args(OptionParser(option_list = option_list))

samples <- trimws(strsplit(opt$samples, ",")[[1]])
samples <- samples[nzchar(samples)]
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

abu <- list()
cnt <- list()
for (s in samples) {
  fa <- file.path(opt$abundance_dir, paste0(s, ".taxon_abundance.tsv"))
  fc <- file.path(opt$counts_dir, paste0(s, ".filtered_counts.tsv"))
  if (!file.exists(fa)) fa <- file.path(opt$abundance_dir, s, paste0(s, ".taxon_abundance.tsv"))
  if (!file.exists(fc)) fc <- file.path(opt$counts_dir, s, paste0(s, ".filtered_counts.tsv"))
  if (file.exists(fa)) {
    x <- fread(fa)
    if (nrow(x) > 0) { x[, sample := s]; abu[[s]] <- x }
  }
  if (file.exists(fc)) {
    y <- fread(fc)
    if (nrow(y) > 0) { y[, sample := s]; cnt[[s]] <- y }
  }
}

if (length(abu) == 0) {
  fwrite(data.table(taxid = character()), file.path(opt$outdir, "taxon_umi_matrix.tsv"),
         sep = "\t", quote = FALSE)
  fwrite(data.table(sample = character(), barcode = character(), taxid = character(),
                    counts = integer()),
         file.path(opt$outdir, "cell_taxon_umi_long.tsv"), sep = "\t", quote = FALSE)
  message("[11] done: no taxa")
  quit(save = "no", status = 0)
}

a <- rbindlist(abu, fill = TRUE)
a[, taxid := as.character(taxid)]
mat <- a %>%
  select(taxid, sample, umi) %>%
  pivot_wider(names_from = sample, values_from = umi, values_fill = 0)
for (s in samples) if (!s %in% names(mat)) mat[[s]] <- 0
mat <- as.data.table(mat)
setcolorder(mat, c("taxid", samples))
fwrite(mat, file.path(opt$outdir, "taxon_umi_matrix.unfiltered.tsv"),
       sep = "\t", quote = FALSE)

if (opt$no_virus_filter) {
  fwrite(mat, file.path(opt$outdir, "taxon_umi_matrix.tsv"), sep = "\t", quote = FALSE)
  keep_taxa <- mat$taxid
  if (length(cnt) > 0) {
    c1 <- rbindlist(cnt, fill = TRUE)
    c1[, taxid := as.character(taxid)]
    keep <- intersect(c("sample", "barcode", "taxid", "counts", "rank", "genus", "species"),
                      names(c1))
    fwrite(c1[taxid %in% keep_taxa, ..keep],
           file.path(opt$outdir, "cell_taxon_umi_long.tsv"), sep = "\t", quote = FALSE)
  } else {
    fwrite(data.table(sample = character(), barcode = character(), taxid = character(),
                      counts = integer()),
           file.path(opt$outdir, "cell_taxon_umi_long.tsv"), sep = "\t", quote = FALSE)
  }
  message("[11] done (viral filter disabled)")
  quit(save = "no", status = 0)
}

if (is.null(opt$lineage) || !file.exists(opt$lineage)) {
  stop("--lineage is required unless --no-virus-filter is used")
}
lineage <- fread(opt$lineage)
required <- c("taxid", "kindom", "genus", "species")
if (!all(required %in% names(lineage))) {
  stop("lineage must contain: ", paste(required, collapse = ", "))
}
lineage[, taxid := as.character(taxid)]
lineage <- unique(lineage[, .(taxid, kingdom = kindom,
                              genus = as.character(genus),
                              species = as.character(species))])
lineage[, species_clean := ifelse(is.na(species) | species == "" | species == "NA",
                                  genus, species)]
lineage[, species_clean := sub("^[sg]__", "", species_clean)]
lineage[, genus_clean := sub("^[sg]__", "", genus)]

ann <- merge(data.table(taxid = mat$taxid), lineage, by = "taxid", all.x = TRUE)
ann[is.na(kingdom), kingdom := "NA"]
ann[is.na(species_clean), species_clean := ""]

excluded_names <- c("Casadabanvirus D3112", "Hepatitis E virus",
                    "Kayfunavirus LM33P1", "Metaavulavirus yucaipaense",
                    "Orthonairovirus khani", "Webervirus KLPPOU149",
                    "Melbournevirus")

read_ncbi_host <- function(path) {
  x <- fread(path, sep = "\t", quote = "")
  req <- c("Organism_Name", "Species", "Genus", "Host")
  if (!all(req %in% names(x))) stop("NCBI host metadata must contain: ",
                                    paste(req, collapse = ", "))
  x <- unique(x[, ..req])
  x <- x[!Organism_Name %in% excluded_names]
  long <- rbind(
    x[!is.na(Organism_Name) & nzchar(Organism_Name), .(name = Organism_Name, host = Host)],
    x[!is.na(Species) & nzchar(Species), .(name = Species, host = Host)],
    x[!is.na(Genus) & nzchar(Genus), .(name = Genus, host = Host)]
  )
  long <- long[!is.na(name) & nzchar(name) & !is.na(host) & nzchar(host)]
  long[, .(host = paste(unique(host), collapse = ";")), by = name]
}

ncbi <- NULL
if (!is.null(opt$ncbi_virus_host) && nzchar(opt$ncbi_virus_host)) {
  if (!file.exists(opt$ncbi_virus_host)) {
    stop("NCBI Virus host metadata not found: ", opt$ncbi_virus_host)
  }
  ncbi <- read_ncbi_host(opt$ncbi_virus_host)
}
if (!is.null(opt$ncbi_virus_assembly_host) && nzchar(opt$ncbi_virus_assembly_host)) {
  if (!file.exists(opt$ncbi_virus_assembly_host)) {
    stop("NCBI Virus assembly host metadata not found: ", opt$ncbi_virus_assembly_host)
  }
  assembly <- fread(opt$ncbi_virus_assembly_host)
  req <- c("Organism_Name", "Species", "Genus", "Host")
  if (all(req %in% names(assembly))) {
    tmp <- tempfile(fileext = ".tsv")
    fwrite(unique(assembly[, ..req]), tmp, sep = "\t", quote = FALSE)
    assembly_map <- read_ncbi_host(tmp)
    unlink(tmp)
    ncbi <- if (is.null(ncbi)) assembly_map else {
      rbind(ncbi, assembly_map)[, .(host = paste(unique(strsplit(paste(host, collapse = ";"), ";")[[1]]),
                                                collapse = ";")), by = name]
    }
  }
}

hvd_names <- character()
if (!is.null(opt$hvd) && nzchar(opt$hvd)) {
  if (!file.exists(opt$hvd)) stop("HVD table not found: ", opt$hvd)
  hvd <- tryCatch(fread(opt$hvd, sep = ","), error = function(e) NULL)
  if (is.null(hvd) || ncol(hvd) == 1) hvd <- fread(opt$hvd, sep = "\t")
  if ("VirusName" %in% names(hvd)) hvd_names <- unique(trimws(as.character(hvd$VirusName)))
}

hvidb_taxids <- character()
if (!is.null(opt$hvidb) && nzchar(opt$hvidb)) {
  if (!file.exists(opt$hvidb)) stop("HVIDB table not found: ", opt$hvidb)
  hvidb <- fread(opt$hvidb, sep = "\t", quote = "")
  if ("Organism_Interactor_virus" %in% names(hvidb)) {
    x <- as.character(hvidb$Organism_Interactor_virus)
    hvidb_taxids <- unique(trimws(unlist(strsplit(x, "[,;[:space:]]+"))))
    hvidb_taxids <- hvidb_taxids[grepl("^[0-9]+$", hvidb_taxids)]
  }
}

ann[, host := "unknown"]
if (!is.null(ncbi)) {
  ann[ncbi, on = .(species_clean = name), host := i.host]
}
ann[, host := ifelse(is.na(host) | host == "", "unknown", host)]
ann[, ncbihost_human := ifelse(grepl("Homo sapiens", host), "YES", "NO")]
ann[, humanvirusbase := ifelse(species_clean %in% hvd_names, "Human", "NO")]
ann[, HVIDB_human := ifelse(taxid %in% hvidb_taxids, "YES", "NO")]
ann[, virus_score := as.integer(ncbihost_human == "YES") +
                    as.integer(humanvirusbase == "Human") +
                    as.integer(HVIDB_human == "YES")]

ann[, kept := kingdom == "k__Bacteria" |
                (kingdom == "k__Viruses" & virus_score >= opt$virus_min_score) |
                (isTRUE(opt$keep_other_kingdoms) &
                   !(kingdom %in% c("k__Bacteria", "k__Viruses")))]
ann[, reason := fifelse(kept, "kept",
                        fifelse(kingdom == "k__Viruses",
                                "virus_low_host_evidence",
                                "other_kingdom"))]
ann[is.na(kingdom) | kingdom == "NA", reason := "missing_lineage"]
ann[is.na(kept), kept := FALSE]

keep_taxa <- ann[kept == TRUE, taxid]
fwrite(ann[, .(taxid, kingdom, genus, species, species_clean, host,
               ncbihost_human, humanvirusbase, HVIDB_human, virus_score, kept, reason)],
       file.path(opt$outdir, "taxon_annotation.tsv"), sep = "\t", quote = FALSE)

mat_keep <- mat[taxid %in% keep_taxa]
fwrite(mat_keep, file.path(opt$outdir, "taxon_umi_matrix.tsv"),
       sep = "\t", quote = FALSE)

if (length(cnt) > 0) {
  c1 <- rbindlist(cnt, fill = TRUE)
  c1[, taxid := as.character(taxid)]
  keep <- intersect(c("sample", "barcode", "taxid", "counts", "rank", "genus", "species"),
                    names(c1))
  fwrite(c1[taxid %in% keep_taxa, ..keep],
         file.path(opt$outdir, "cell_taxon_umi_long.tsv"), sep = "\t", quote = FALSE)
} else {
  fwrite(data.table(sample = character(), barcode = character(), taxid = character(),
                    counts = integer()),
         file.path(opt$outdir, "cell_taxon_umi_long.tsv"), sep = "\t", quote = FALSE)
}

summary_dt <- ann[, .(n_taxa = .N, n_kept = sum(kept)), by = .(kingdom, reason)]
fwrite(summary_dt, file.path(opt$outdir, "virus_filter_summary.tsv"),
       sep = "\t", quote = FALSE)

message("[11] done: ", length(keep_taxa), " retained taxa")
