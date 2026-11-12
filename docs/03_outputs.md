# Output reference

All paths are relative to `OUTDIR` (default: `results/`).

## Main deliverables

### `11_summary/taxon_umi_matrix.tsv`

Taxon × sample matrix. First column is `taxid`; each remaining column is one
sample; values are summed UMI counts for the retained taxa.

```text
taxid	sample1	sample2	sample3
1000	40	80	120
1002	20	40	60
```

### `11_summary/cell_taxon_umi_long.tsv`

One row per sample × barcode × taxon:

```text
sample	barcode	taxid	counts	rank	genus	species
sample1	AAACCCGGGTTTAAAA	1000	5	g	ToygenusA	
```

This is the table to join to Seurat `meta.data` or to a cell-level metadata
table. `barcode` is the 16-bp barcode without the `-1` suffix; step 10 also
checks the Cell Ranger barcode file when `CELLRANGER_FILTER=1`.

## Per-sample outputs

### `10_quantify/<sample>/<sample>.qc.tsv`

```text
sample	microbial_umi	taxa_detected	barcodes_with_microbes
```

Useful as a first sanity check: if `taxa_detected` is zero, inspect
`taxon_metrics.tsv` and the earlier steps.

### `10_quantify/<sample>/<sample>.taxon_metrics.tsv`

One row per candidate taxon with the values used by the final filter:

- `rpmm`: taxon reads per million microbial reads;
- `CLrpmm`: cell-line 90th-percentile threshold;
- `r`: barcode-level k-mer correlation;
- `r1`, `r2`, `r3`: sample-level correlations;
- `p`, `p1`, `p2`, `p3`: corresponding p-values;
- `rank`, `name`, `mreads`, `genus`, `species`;
- `decontam_p`;
- `black` and `pass`.

Look here first when a taxon unexpectedly disappears.

### `10_quantify/<sample>/<sample>.filtered_counts.tsv`

Final barcode × taxon counts after all filters. Columns include barcode,
taxid, rank, lineage, counts, and the filter metrics.

### `10_quantify/<sample>/<sample>.taxon_abundance.tsv`

Taxon-level summary for the sample:

```text
taxid	rank	name	genus	species	umi	...
```

### `10_quantify/<sample>/<sample>.counts.tsv`

Unfiltered barcode × taxon counts after UMI collapsing, before the final
taxon-level filter. Useful for debugging.

## Denoising outputs

### `06_denoise/<sample>/<sample>.barcode.kmer.hits.tsv`

Per-taxon barcode-level correlation:

```text
taxid	r	p	name
```

`r` is the Spearman correlation between total and unique k-mers across
barcodes. It is `NA` when there are too few barcodes.

### `06_denoise/<sample>/<sample>.kraken.report.rpmm.tsv`

Taxon-level report with minimizer/read counts and RPMM. Step 07 uses this file
across samples.

### `07_sample_denoise/sample_denoising.tsv`

```text
taxid	r1	r2	r3	p1	p2	p3	rank	name	nn
```

`nn` is the number of samples in which the taxon was observed. With fewer than
3 samples, the sample-level correlations are under-powered.

## Kraken2/BWA intermediates

| Path | Meaning |
|---|---|
| `02_kraken/<sample>/<sample>.kraken.output.txt` | per-read Kraken2 assignment |
| `02_kraken/<sample>/<sample>.kraken.report.txt` | full report with minimizer columns |
| `02_kraken/<sample>/<sample>.kraken.report.std.txt` | standard 6-column report |
| `02_kraken/<sample>/<sample>.kraken.report.mpa.txt` | MPA lineage table |
| `03_microbial/<sample>/<sample>_1.fa` | microbial R1 FASTA (barcode+UMI) |
| `03_microbial/<sample>/<sample>_2.fa` | microbial R2 FASTA |
| `04_bwa/<sample>/<sample>_bwa_pass_reads.tsv` | reads passing BWA validation |
| `04_bwa/<sample>/<sample>_bwa_black_species_id.tsv` | species taxids that failed BWA validation |
| `05_host_filter/<sample>/<sample>_clean_1.fq` | host-filtered R1 |
| `05_host_filter/<sample>/<sample>_clean_2.fq` | host-filtered R2 |

## Logs

`logs/<sample>.<step>.log` contains the combined stdout/stderr for a sample
step. `logs/07_sample_denoising.log`, `logs/09_decontam.log`, and
`logs/11_summary.log` cover project-level steps.
