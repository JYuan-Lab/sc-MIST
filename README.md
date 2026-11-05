
# sc-MIST
![sc-MIST workflow](assets/sc-MIST_pipeline.png)

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
4. False-positive filtering
        |
        v
5. Contamination filtering
        |
        v
6. Microbial quantification

```
---


## Reference resources

sc-MIST requires the following reference resources. Large reference databases are not included in this repository; users should download or build them from the sources below.

1. Human reference for host alignment

For 10x Genomics scRNA-seq data, host alignment was performed using the Cell Ranger GRCh38-2020-A reference (GENCODE v32 / Ensembl 98).


2. Custom Kraken2 database

The sc-MIST Kraken2 database consists of:

Human: T2T-CHM13v2.0
Bacteria: Kraken2 bacterial reference genomes, restricted to taxa with prior evidence of human association
Viruses: human-infecting viral genomes from NCBI Virus / RefSeq
Taxonomy: NCBI taxonomy

The bacterial reference was curated by comparing Kraken2 bacterial taxa with human-associated microbial resources including HMP, MGnify, mBodyMap, and published literature (PMID: 32467386).

The overlap and scoring results used for bacterial reference curation are provided in [`resources/kraken_overlap_databases.xlsx`](resources/kraken_overlap_databases.xlsx).

Matches were scored by taxonomic level:

species = 3
genus   = 2
phylum  = 1

Species with a cumulative score >= 6 were retained, resulting in 3,808 bacterial species in the curated reference.

The custom database can then be constructed using standard Kraken2 database-building procedures.

3. Taxon-specific BWA reference library

A separate BWA reference library is used to validate Kraken2 assignments.

Species-level assignment: reference sequences from the corresponding species are indexed directly.
Genus-level assignment: sequences from retained species belonging to the same genus are combined and indexed together.

Bacterial and viral taxa included in the sc-MIST microbial reference should have corresponding BWA indices.

4. HLA reference

Residual reads derived from highly polymorphic HLA/MHC regions are filtered using genomic HLA sequences from the IPD-IMGT/HLA Database.

5. Human rRNA reference

Residual human ribosomal RNA reads are filtered using an NCBI-derived reference containing:

18S rRNA
28S rRNA
mitochondrial rRNA


---

# 1. Human-unmapped read extraction

scRNA-seq reads are first aligned to the human reference genome to remove host-derived reads.

* **10x Genomics:** Cell Ranger
* **non-10x droplet-based scRNA-seq:** zUMI

> **Note:** For 10x Genomics datasets, we recommend running Cell Ranger before starting the downstream sc-MIST workflow. By default, the provided Snakemake workflow assumes that Cell Ranger processing has already been completed and therefore does not run Cell Ranger automatically.

Only reads that remain unmapped to the human reference are retained for downstream microbial analysis.

### Input

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
human reference genome
```

For 10x Genomics data, the downstream workflow can also start from a Cell Ranger output BAM file:

```text
possorted_genome_bam.bam
```

### Output

```text
sample.unmapped.R1.fastq.gz
sample.unmapped.R2.fastq.gz
```

For droplet-based datasets, cell barcode and UMI information are preserved so that microbial reads can later be reassigned to individual cells.

### Human reference

For 10x Genomics data, our implementation used the Cell Ranger GRCh38 reference:

```text
refdata-gex-GRCh38-2020-A
```


---

# 2. Initial microbial classification

Human-unmapped reads are classified using **Kraken2/Kraken2Uniq** against the custom sc-MIST reference database described in [Reference resources](#reference-resources).

Kraken2 assigns candidate taxonomic labels based on k-mer matches and the lowest common ancestor (LCA) framework.

The custom database contains human reference sequences, curated human-associated bacterial taxa, and human-infecting viral sequences.

### Input

```text
human-unmapped FASTQ
custom sc-MIST Kraken2 database
```

### Output

```text
sample.kraken.output.txt
sample.kraken.report.txt
sample.kraken.report.mpa.txt
```

The Kraken output provides read-level taxonomic assignments and taxon-level summaries that are used for downstream sequence validation.

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

---

# 4. False-positive filtering

sc-MIST applies multiple denoising criteria to retain taxa with internally consistent sequence support.


## 4.1 Barcode-level correlation

For droplet-based datasets, total k-mers and unique k-mers are summarized across cell barcodes for each taxon.

A Spearman correlation is calculated between:

```text
total k-mers
vs
unique k-mers
```

Taxa are retained when:

```text
r > 0
```

---

## 4.2 Sample-level correlation

Within each cohort, sample-level correlation is performed when:

```text
number of samples >= 2
```

For each taxon, three Spearman correlations are calculated:

```text
r1:
unique k-mers
vs
total k-mers

r2:
total reads
vs
total k-mers

r3:
unique k-mers
vs
total reads
```

Taxa are retained when:

```text
r1 > 0
r2 > 0
r3 > 0
```

---



# 5. Decontamination

## 5.1 Cell-line background filtering

For each taxon, microbial abundance is normalized as RPMM:

```text
RPMM =
taxon reads / total microbial reads x 1,000,000
```

Sample RPMM is compared with the corresponding taxon-specific background distribution from cell-line datasets.

Current threshold:

```text
Taxid RPMM > cell-line 90th percentile RPMM
```


---

## 5.2 Decontam

The R package **decontam** is used with the frequency-based method to identify potential contaminants.

Parameters:

```text
method: frequency
threshold: [TO ADD]
concentration variable: [TO ADD]
```

---

## 5.3 Contaminant blacklists

Additional filtering is performed using:

1. a literature-curated laboratory contaminant blacklist：
   ([`resources/black_list.txt`](resources/black_list.txt));
2. an MGnify-derived environmental blacklist.
   ([`resources/mgnify_black.txt`](resources/mgnify_black.txt)).

For the MGnify blacklist, taxa reported in more than four environment-associated studies can be flagged.

---

# 6. Microbial quantification

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

---

# Repository structure

```text
sc-MIST/
├── README.md
├── LICENSE
├── .gitignore
├── scMIST_pipeline.smk
├── config.example.yaml
├── environment.yaml
│
├── scripts/
│   └── barcode_denoising_v5.r
│   └── bwa_filter_v2.sh
│   └── bwa_mhc_rrna.r
│   └── ...
│
└── resources/
    ├── black_list.txt
    ├── kraken_overlap_databases.xlsx
    └── mgnify_black.txt


```

---

# Configuration

Example:

```yaml
config_path: /users/project_name/config/config_example.yaml
Output_dir: /users/project_name/output
scMIST: /users/project_name/software/scMIST
fq_path: /users/project_name/dataset/cohort
ncbi_blast: /users/project_name/software/ncbi-blast/bin
Kraken2Uniq_path: /users/project_name/software/kraken2/kraken2
kraken_library_path: /users/project_name/reference/kraken_database/database_ncbi_T2T_asshuman_score_6
umi: 10  #v2 10 v3 12
cb: 16  #v2 16 v3 16
R1: 1
R2: 2
sample: sample1, sample2 ,sample3, sample4
```


---

# Running sc-MIST


Run:

```bash
conda env create -f environment.yaml -n sc-MIST
conda activate sc-MIST

snakemake -s sc-MIST_pipeline.smk \
    --cores 8 \
    --configfile config.example.yaml

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
