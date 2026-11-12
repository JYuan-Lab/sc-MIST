# Troubleshooting

## Input and environment

### `bamtofastq: command not found`

The pipeline uses 10x `bamtofastq` to convert the Cell Ranger BAM. Install
Cell Ranger and add its `bamtofastq` executable to `PATH`, or set:

```bash
BAMTOFASTQ=/absolute/path/to/bamtofastq
```

The demo uses `demo/utils/bamtofastq_samtools.sh`; do not use that wrapper for
real Cell Ranger BAMs unless you have checked that it preserves your read
structure.

### `Rscript: command not found`

Activate the conda environment:

```bash
conda activate scmist
which Rscript
```

### R package missing

Check:

```bash
Rscript -e 'for (p in c("optparse","dplyr","tidyr","stringr","data.table","readr","Matrix","decontam")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'
```

Install a missing Bioconductor package with `BiocManager::install()`. Recreate
the environment from `environment.yaml` if the environment is badly broken.

## Sample and threshold problems

### Fewer than 3 samples

Step 07 needs at least 3 samples. The runner prints a warning but continues;
the strict filter will be unreliable. Add samples rather than lowering the
threshold.

### No taxa in `taxon_umi_matrix.tsv`

Open the per-sample metrics:

```text
results/10_quantify/<sample>/<sample>.taxon_metrics.tsv
```

Check each condition in `pass`. The most common causes are:

- `r`, `r1`, `r2`, or `r3` is `NA` or below `MIN_CORR`;
- `rpmm <= CLrpmm`;
- the taxon is in a blacklist;
- `decontam_p <= 0.1`;
- the taxon is a species in the BWA blacklist.

For a first pass, try a small read-depth dataset or a less stringent
`MIN_CORR`, but report the change.

### All taxa are `NA` in step 06

Step 06 needs multiple barcodes per taxon. A library with only one or two
barcodes cannot estimate a barcode-level correlation. This is expected for
tiny test data; use more cells/barcodes for real data.

### `decontam` not installed

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("decontam")
```

`--skip-decontam` writes `NA` p-values and is only for demos/tests.

### Seurat v5 errors

The pipeline scripts do not use Seurat. If you use Seurat downstream, use v4.4.0.
Check and reinstall:

```r
packageVersion("Seurat")
remotes::install_version("Seurat", version = "4.4.0",
                         repos = "https://cloud.r-project.org")
```

## Reference problems

### `Kraken2 database errors`

The database directory must contain:

```text
hash.k2d
opts.k2d
taxo.k2d
```

Also confirm that the database k-mer length matches the one used to build it.
The released database used `k=35`.

### `Missing BWA reference: .../<taxid>.fa`

The BWA library does not contain an index for a taxid present in the lineage
matrix. Build the missing species/genus FASTA and run `bwa index`. The helper
`scripts/utils/build_bwa_library.sh` can rebuild the whole library from
`<taxid>.fa` files.

### `BWA index not found for host_residual.fa`

```bash
bwa index "$HOST_RESIDUAL_FASTA"
```

### Lineage matrix errors

The required columns are:

```text
taxid	kindom	phylum	genus	species	g_taxid	s_taxid
```

The `kindom` spelling is intentional. `NA` is allowed for missing ranks.

### Cell-line table errors

`08_filter_cellline.R` accepts either `q90` or
`threshold_0.9_exclude0`. The file must contain a `taxid` column. If the table
is missing, use a table with `q90 = 0` only to test the mechanics; this
disables the cell-line filter.

## Resuming and overwriting

`run_pipeline.sh` does not delete old output. If you change a threshold, rerun
the affected step and all downstream steps. For example, changing step 07
requires rerunning 07, 08, 09, 10, and 11.

```bash
bash scripts/run_pipeline.sh --config config/config.sh --start 7 --end 11
```

Keep a copy of a successful run before experimenting with thresholds.
