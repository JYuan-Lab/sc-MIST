# Expected demo output

These files are reference outputs from the synthetic demo. They let you check
that the pipeline ran correctly without comparing every intermediate file.

| File | What to check |
|---|---|
| `taxon_umi_matrix.tsv` | four retained taxa: 1000, 1002, 2001, 3001 |
| `qc_summary.tsv` | 4 taxa detected in each of the 3 samples |
| `cell_taxon_umi_long.tsv` | cell-level counts for the retained taxa |
| `sha256sums.txt` | checksums of the three tables above |

After a run, compare the observed tables with these files:

```bash
bash demo/check_results.sh
```

Expected matrix:

```text
taxid	sample1	sample2	sample3
1000	40	80	120
1002	20	40	60
2001	20	40	60
3001	20	40	60
```

Expected QC:

```text
sample	microbial_umi	taxa_detected	barcodes_with_microbes
sample1	100	4	8
sample2	200	4	8
sample3	300	4	8
```

Exact counts can shift by a read or two with different Kraken2/BWA versions.
The retained taxid set and the non-empty three-sample matrix are the important
checks. The demo runs DECONTAM when the package is available; if it is missing,
step 09 writes `NA` p-values and the same expected tables are produced.
