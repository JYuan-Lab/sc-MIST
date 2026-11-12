# Step-by-step reference

The full commands and runnable examples are in the main `README.md`. This page
is a compact map of the 11 scripts, their inputs, and their outputs.

| Step | Script | Input | Key output |
|---|---|---|---|
| 01 | `01_extract_unmapped.sh` | Cell Ranger `possorted_genome_bam.bam` | unmapped BAM; R1/R2 FASTQ; R1 barcode+UMI FASTQ |
| 02 | `02_kraken_classify.sh` | step 01 FASTQ; Kraken2 DB | Kraken output/report/MPA; classified FASTQ |
| 03 | `03_extract_microbial_reads.R` | step 02 reports/FASTQ | microbial FASTA; microbiome output |
| 04 | `04_validate_taxa_bwa.R` | step 03 output; lineage; curated species; BWA library | BWA-pass reads; BWA blacklist |
| 05 | `05_filter_host_and_refresh.sh` | step 04 FASTQ; HLA/rRNA FASTA; Kraken2 DB | clean FASTQ; refreshed Kraken reports; microbial FASTA |
| 06 | `06_barcode_denoising.R` | step 05 output; step 02 report | barcode k-mer hits; RPMM table |
| 07 | `07_sample_denoising.R` | all step 06 RPMM tables | `sample_denoising.tsv` |
| 08 | `08_filter_cellline.R` | step 06 hits/RPMM; step 07 table; step 05 report; q90 table | cell-line hit tables |
| 09 | `09_decontam.R` | all step 08 hits; step 02 reports | `decontam.tsv` |
| 10 | `10_quantify_filter.R` | step 05 FASTA; step 08 hits; step 09 table; blacklists | counts; filtered counts; abundance; QC |
| 11 | `11_summarize.R` | all step 10 tables; lineage; viral evidence | taxon × sample and cell × taxon matrices |

`05_filter_host_residuals.R` is an internal helper called by
`05_filter_host_and_refresh.sh`; it is not a separate numbered step.

## Execution order in `run_pipeline.sh`

The runner does **not** execute 01→11 inside one per-sample loop. Step 07 needs
all samples, and step 09 needs all step 08 outputs. The actual order is:

```text
for each sample: 01 → 06
07 (all samples)
for each sample: 08
09 (all samples)
for each sample: 10
11
```

This is why `--start`/`--end` should be used with care. If you rerun only step
11, the step 10 tables must already exist. If you rerun step 09, all step 08
hits must already exist.

## Minimal per-step invocation

Each script has a `--help` flag:

```bash
bash scripts/01_extract_unmapped.sh --help
bash scripts/02_kraken_classify.sh --help
Rscript scripts/03_extract_microbial_reads.R --help
Rscript scripts/04_validate_taxa_bwa.R --help
bash scripts/05_filter_host_and_refresh.sh --help
Rscript scripts/05_filter_host_residuals.R --help
Rscript scripts/06_barcode_denoising.R --help
Rscript scripts/07_sample_denoising.R --help
Rscript scripts/08_filter_cellline.R --help
Rscript scripts/09_decontam.R --help
Rscript scripts/10_quantify_filter.R --help
Rscript scripts/11_summarize.R --help
```

## Resuming

```bash
# Final quantification and project summary only
bash scripts/run_pipeline.sh --config config/config.sh --start 10 --end 11

# Recompute sample-level, cell-line, and contamination filters for an existing run
bash scripts/run_pipeline.sh --config config/config.sh --start 7 --end 9
```

The runner does not delete previous output. Files are overwritten in place, so
keep a copy of a successful run before changing thresholds.
