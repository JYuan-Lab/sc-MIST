# Quick start

This is the shortest path from a fresh checkout to a demo result and then to a
real project.

## 1. Install

```bash
conda env create -f environment.yaml
conda activate scmist
```

Install Cell Ranger separately if you need to produce `outs/` from FASTQ. The
pipeline starts from `outs/possorted_genome_bam.bam`, so Cell Ranger does not
need to be on the same machine once the `outs/` directories exist.

Check R:

```bash
Rscript -e 'for (p in c("optparse","dplyr","tidyr","stringr","data.table","readr","Matrix","decontam")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'
```

The manuscript used Kraken2 2.1.3, BWA 0.7.18, samtools 1.21, SeqKit 2.6.1,
taxonkit 0.16.0, and R 4.4. `environment.yaml` uses portable version bounds;
pin exact versions on Linux only if you need a byte-for-byte reproduction.

## 2. Run the demo

```bash
bash demo/run_demo.sh
cat demo/results/11_summary/taxon_umi_matrix.tsv
bash demo/check_results.sh
```

Expected retained taxa: `1000`, `1002`, `2001`, `3001` across three samples.
See `demo/expected/` for reference tables.

## 3. Prepare a real project

```bash
cp config/config.example.sh config/config.sh
$EDITOR config/config.sh
cp config/samples.example.tsv config/samples.tsv
$EDITOR config/samples.tsv
```

Fill in:

- `KRAKEN_DB`
- `BWA_LIBRARY`
- `TAXON_LINEAGE`
- `HOST_RESIDUAL_FASTA`
- `CELLLINE_Q90`
- one Cell Ranger `outs/` path per sample (at least 3 samples)

Build the large references first (see `docs/01_reference_build.md`), then run
the checker:

```bash
bash scripts/utils/check_references.sh --config config/config.sh
```

## 4. Run

```bash
bash scripts/run_pipeline.sh --config config/config.sh
```

Start from a later step after fixing an input:

```bash
bash scripts/run_pipeline.sh --config config/config.sh --start 9 --end 11
```

## 5. Inspect

```text
results/10_quantify/<sample>/<sample>.taxon_metrics.tsv   # why taxa passed/failed
results/11_summary/taxon_umi_matrix.tsv                  # taxon x sample
results/11_summary/cell_taxon_umi_long.tsv               # cell-level
```

For reference construction, see `docs/01_reference_build.md`. For every step's
inputs, commands, and outputs, see the main `README.md` and
`docs/02_step_by_step.md`.
