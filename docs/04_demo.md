# Demo walkthrough

The demo is synthetic. It is designed to exercise all 11 steps with a tiny
reference and three samples, not to make biological claims.

## What is generated

`demo/make_demo.py` creates:

```text
demo/input/
├── sample1/outs/
│   ├── possorted_genome_bam.bam
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
├── sample2/outs/...
└── sample3/outs/...
```

The toy reference source files are committed. Binary Kraken2 and BWA indices
are rebuilt by `demo/build_reference.sh` on your machine, because those formats
are platform-dependent. The toy reference contains:

```text
taxid 1000  ToygenusA       (genus)
taxid 1001  ToyspeciesA     (species)
taxid 1002  ToyspeciesA2    (species)
taxid 2000  ToygenusB       (genus)
taxid 2001  ToyspeciesB     (species)
taxid 3000  Toyvirusgenus   (genus)
taxid 3001  Toyvirus        (species)
taxid 9606  Homo sapiens    (host)
```

Species 1001 and 1002 share a 250-bp block. Reads from that block are assigned
to genus 1000 by the LCA rule; reads from the unique regions are assigned to
the species. This demonstrates why both genus- and species-level references
are needed.

## Run it

```bash
conda activate scmist
bash demo/run_demo.sh
```

The script does three things:

1. `python3 demo/make_demo.py` writes the toy Cell Ranger output and source
   reference FASTA/tables.
2. `bash demo/build_reference.sh` builds the toy Kraken2 database and BWA
   indices.
3. `bash scripts/run_pipeline.sh --config demo/config.demo.sh --start 1 --end 11`
   runs the pipeline. Step 09 uses the real DECONTAM frequency model when the
   package is installed and falls back to `NA` p-values with a warning
   otherwise.

Then compare the observed tables with the reference tables:

```bash
bash demo/check_results.sh
```

The demo uses `demo/utils/bamtofastq_samtools.sh`, a small `samtools`-based
replacement for the 10x `bamtofastq` executable. This is only for portability;
real Cell Ranger BAMs should use the official 10x tool.

## Expected result

`demo/results/11_summary/taxon_umi_matrix.tsv`:

```text
taxid	sample1	sample2	sample3
1000	40	80	120
1002	20	40	60
2001	20	40	60
3001	20	40	60
```

With only three samples, DECONTAM prints a warning that some batches have very
few samples. This is expected for the demo; use more samples in a real project.
Four taxa pass the strict filters. `1001` is removed by the host-residual filter
because the toy HLA reference intentionally contains part of its sequence;
`2000` and `3000` are genus taxids without barcode-level hits. This is expected
and shows that the filters are active.

Observed values for `sample1` (project-level steps use all three samples):

| Step | Main output | Demo observation |
|---:|---|---|
| 01 | unmapped BAM + R1/R2 + barcode+UMI R1 | 103 read pairs |
| 02 | Kraken output/report/MPA + classified FASTQ | 103 pairs classified |
| 03 | microbial FASTA + microbiome output | 101 R1 records; 11 microbial taxids |
| 04 | BWA-pass reads + BWA blacklist | 7 validated taxa; 182 pass records |
| 05 | clean FASTQ + refreshed Kraken output | 80 clean pairs; 10 refreshed taxids |
| 06 | barcode k-mer hits + RPMM table | 4 taxa with barcode hits |
| 07 | `sample_denoising.tsv` | 6 taxa across 3 samples |
| 08 | cell-line hit tables | 6 candidate taxa per sample |
| 09 | `decontam.tsv` | 6 taxa; none flagged (p ~ 1) |
| 10 | counts/filtered counts/abundance/QC | 4 retained; 100 UMIs; 8 barcodes |
| 11 | final matrices and annotations | 4 taxa × 3 samples |

Exact values may vary slightly with Kraken2/BWA versions. The key expected
properties are:

- all three samples are represented;
- the matrix is non-empty;
- the retained taxid set is `1000, 1002, 2001, 3001`;
- `demo/results/10_quantify/<sample>/<sample>.qc.tsv` reports
  `taxa_detected = 4`.

## Rebuild only the references

```bash
bash demo/build_reference.sh
```

## Re-run only the pipeline

```bash
bash scripts/run_pipeline.sh \
  --config demo/config.demo.sh \
  --start 1 --end 11
```

## Clean up demo results

The demo writes to `demo/results/`; it is ignored by `.gitignore`. To start
fresh, remove that directory with your usual file manager or `find` command.
`demo/input/` and `demo/reference/` can be regenerated with
`python3 demo/make_demo.py` and `bash demo/build_reference.sh`.
