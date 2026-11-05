![sc-MIST workflow](assets/scMIST_workflow.png)

# sc-MIST

**sc-MIST** (**single-cell Microbial Identification with Stringent Testing**) is a workflow for recovering high-confidence microbial signals from host-dominated single-cell RNA-seq (scRNA-seq) data while reducing spurious taxonomic assignments and contamination.

sc-MIST was developed for the study:

> **A Pan-Cancer Cell-Associated Microbiome Atlas Reveals Inter- and Intra-Tumor Heterogeneity**
> *Under Review*

---

## Overview

sc-MIST combines host-read depletion, microbial classification, sequence-level validation, statistical denoising, contamination filtering, and barcode-aware microbial quantification.

To reduce false-positive assignments, the microbial reference is restricted to taxa with prior evidence of human association. The curated reference includes human-associated bacteria from MGnify, mBodyMap, the Human Microbiome Project, published literature, and NCBI-annotated human-infecting viruses.

---

# Workflow

```text
Raw scRNA-seq FASTQ
        |
        v
1. Human-unmapped read extraction
        |
        v
2. Initial microbial classification
        |
        v
3. BWA validation and host-residual removal
        |
        v
4. Kraken2 re-classification
        |
        v
5. False-positive filtering
        |
        v
6. Contamination filtering
        |
        v
7. Microbial quantification
        |
        v
High-confidence Cell x Microbe matrix
```

---

# 1. Human-unmapped read extraction

scRNA-seq reads are first aligned to the human reference genome.

* **10x Genomics:** Cell Ranger
* **non-10x droplet-based scRNA-seq:** zUMI

Only reads remaining unmapped to the human reference are retained.

### Input

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
human reference genome
```

### Output

```text
sample.unmapped.R1.fastq.gz
sample.unmapped.R2.fastq.gz
```

For droplet-based datasets, cell barcode and UMI information are preserved.

### Human reference

In our implementation, the **T2T-CHM13v2.0** human reference genome was used.

Reference:

[TO ADD: human reference URL / version details]

---

# 2. Initial microbial classification

Human-unmapped reads are classified using **Kraken2/Kraken2Uniq** against a custom microbial reference database.

Each read is assigned a putative taxonomic label using k-mer evidence and the lowest common ancestor framework.

### Custom reference

The database includes:

* human reference sequences;
* human-associated bacterial taxa;
* human-infecting viruses.

Human-associated bacteria were curated from:

* MGnify;
* mBodyMap;
* Human Microbiome Project;
* published literature.

Human-infecting viruses were obtained from NCBI Virus.

### Input

```text
human-unmapped FASTQ
custom Kraken2 database
```

### Output

```text
sample.kraken.output.txt
sample.kraken.report.txt
sample.kraken.report.mpa.txt
```

### Database details

```text
Kraken2 version: [TO ADD]
Kraken2Uniq version: [TO ADD]
k-mer length: 35
NCBI taxonomy version/date: [TO ADD]
viral database date: [TO ADD]
bacterial reference source/version: [TO ADD]
```

---

# 3. BWA validation and host-residual removal

## 3.1 Taxon-specific BWA validation

Reads assigned to microbial species or genera are aligned to the corresponding taxon-specific reference using **BWA-MEM**.

Reads are retained when:

```text
alignment length > 45 bp
sequence identity > 95%
```

### Input

```text
Kraken-positive microbial reads
read-to-taxid assignments
taxon-specific BWA reference library
```

### Output

```text
sample.bwa_pass.R1.fastq.gz
sample.bwa_pass.R2.fastq.gz
```

---

## 3.2 Removal of MHC/HLA and human rRNA reads

BWA-validated reads are further aligned against dedicated human references containing:

* MHC/HLA sequences;
* human rRNA;
* mitochondrial rRNA where applicable.

Reads aligning to these references are removed.

### Input

```text
BWA-validated candidate reads
MHC/HLA reference
human rRNA reference
```

### Output

```text
sample.clean.R1.fastq.gz
sample.clean.R2.fastq.gz
```

Reference details:

```text
MHC/HLA reference: [TO ADD]
human rRNA reference: [TO ADD]
```

---

# 4. Kraken2 re-classification

After BWA validation and host-residual removal, the retained reads are re-classified using Kraken2.

This step regenerates taxonomic assignments for the cleaned read pool and provides the required input for downstream denoising.

### Input

```text
sample.clean.R1.fastq.gz
sample.clean.R2.fastq.gz
custom Kraken2 database
```

### Output

```text
sample.clean.kraken.output.txt
sample.clean.kraken.report.txt
sample.clean.kraken.report.mpa.txt
```

---

# 5. False-positive filtering

sc-MIST applies multiple denoising criteria to retain taxa with internally consistent sequence support.

## 5.1 Read-level k-mer support

A read is retained when:

```text
>= 50% of informative k-mers directly support the assigned taxon
```

---

## 5.2 Barcode-level correlation

For droplet-based datasets, total k-mers and unique k-mers are summarized across cell barcodes for each taxon.

A Spearman correlation is calculated between:

```text
total k-mers
vs
unique k-mers
```

Taxa are retained when:

```text
r > 0.1
```

This step is not applied to non-droplet data.

---

## 5.3 Sample-level correlation

Within each cohort, sample-level correlation is performed when:

```text
number of samples > 2
```

For each taxon, three Spearman correlations are calculated:

```text
r1:
unique minimizers
vs
total minimizers

r2:
total reads
vs
total minimizers

r3:
unique minimizers
vs
total reads
```

Taxa are retained when:

```text
r1 > 0.1
r2 > 0.1
r3 > 0.1
```

---

## 5.4 Genome-distribution assessment

Candidate microbial reads are evaluated for their distribution across the assigned microbial genome.

Broadly distributed reads are considered more consistent with a true microbial signal, whereas highly localized read pile-ups may indicate spurious cross-mapping.

Filtering rule:

```text
[TO ADD: exact genome-distribution metric and threshold]
```

---

# 6. Contamination filtering

## 6.1 Cell-line background filtering

For each taxon, microbial abundance is normalized as RPMM:

```text
RPMM =
taxon reads / total microbial reads x 1,000,000
```

Sample RPMM is compared with the corresponding taxon-specific background distribution from cell-line datasets.

Current threshold:

```text
sample RPMM > cell-line 90th percentile RPMM
```

Reference dataset:

```text
[TO ADD: CCLE/cell-line dataset description and source]
```

---

## 6.2 Decontam

The R package **decontam** is used with the frequency-based method to identify potential contaminants.

Parameters:

```text
method: frequency
threshold: [TO ADD]
concentration variable: [TO ADD]
```

---

## 6.3 Contaminant blacklists

Additional filtering is performed using:

1. a literature-curated laboratory contaminant blacklist;
2. an MGnify-derived environmental blacklist.

For the MGnify blacklist, taxa reported in more than four environment-associated studies can be flagged.

Blacklist resources:

```text
laboratory contaminant reference: [TO ADD]
MGnify filtering procedure: [TO ADD]
```

---

# 7. Microbial quantification

After all filtering steps, retained microbial reads are linked back to individual cells using their original cell barcodes.

For UMI-based scRNA-seq, microbial abundance is quantified using unique UMIs.

Reads are grouped by:

```text
cell barcode
+
microbial taxon
+
UMI
```

PCR duplicates sharing the same barcode, taxon, and UMI are counted once.

### Output

The primary output is a:

```text
Cell x Microbe UMI matrix
```

Example:

| Cell barcode | H_pylori | EBV | HPV16 |
| ------------ | -------: | --: | ----: |
| BC01         |        2 |   0 |     0 |
| BC02         |        0 |   1 |     0 |
| BC03         |        4 |   0 |     2 |

This matrix can be integrated with the corresponding scRNA-seq gene-expression matrix and cell-type annotations.

---

# Input

For 10x Genomics scRNA-seq:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

Typical structure:

```text
R1 = cell barcode + UMI
R2 = transcript sequence
```

For other droplet-based platforms, the appropriate barcode and UMI structure must be specified.

---

# Repository structure

```text
scMIST/
├── README.md
├── LICENSE
├── .gitignore
├── scMIST_pipeline.smk
├── config.example.yaml
├── environment.yaml
│
├── scripts/
│   └── xxx
│   └── xxx
│   └── xxx
│   └── xxx
│   └── xxx
│
└── resources/
    ├── lab_contaminant_blacklist.tsv
    ├── mgnify_environment_blacklist.tsv
    └── human_associated_taxa.tsv


```

---

# Configuration

Example:

```yaml
config_path: /users/project_name/config/config_example.yaml
Output_dir: /users/project_name/output
SAHMI: /users/project_name/software/SAHMI
fq_path: /users/project_name/dataset/cohort
ncbi_blast: /users/project_name/envs/Scasa/bin
Kraken2Uniq_path: /users/project_name/software/kraken2/kraken2
kraken_database_path: /users/project_name/reference/kraken_database/database_ncbi_T2T
umi: 10  #v2 10 v3 12
cb: 16  #v2 16 v3 16
R1: 1
R2: 2
sample: sample1, sample2 ,sample3, sample4
kraken_library_path: /users/project_name/reference/kraken_database/database_ncbi_T2T_asshuman_score_6
```


---

# Running sc-MIST


Run:

```bash
conda env create -f environment.yml -n sc-MIST
conda activate sc-MIST

snakemake -s sc-MIST_pipeline.smk \
    --cores 8 \
    --configfile config.example.yaml \
    --rerun-triggers mtime \
    --rerun-incomplete

echo 'finish'
```


---

# Citation

sc-MIST was developed for:

**A Pan-Cancer Cell-Associated Microbiome Atlas Reveals Inter- and Intra-Tumor Heterogeneity**
*Under Review*

---

# Acknowledgements

sc-MIST builds on ideas from prior work in single-cell microbial signal detection and k-mer-based denoising. We thank the developers of SAHMI and the broader open-source microbial genomics community for providing useful methodological foundations and software resources.

---

# License

This project is distributed under the **MIT License**.



# Contact

Correspondence Email：yuan_jiao@gzlab.ac.cn

---
