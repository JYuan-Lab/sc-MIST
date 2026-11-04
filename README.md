# sc-MIST

**sc-MIST** (**single-cell Microbial Identification with Stringent Testing**) is a workflow for recovering high-confidence microbial signals from host-dominated single-cell RNA-seq (scRNA-seq) data while minimizing spurious taxonomic assignments, host-derived cross-mapping, and experimental contamination.

sc-MIST combines host-read depletion, curated microbial references, k-mer-based taxonomic classification, taxon-specific sequence validation, genome-distribution assessment, statistical denoising, contamination filtering, and barcode-aware microbial quantification.

The workflow was developed for the study:

> **A Pan-Cancer Cell-Associated Microbiome Atlas Reveals Inter- and Intra-Tumor Heterogeneity**
> *Under Review*

---

## Overview

Microbial reads in human scRNA-seq datasets are typically extremely sparse relative to host-derived reads. Their detection is further complicated by several sources of false-positive signal, including:

* incomplete removal of human-derived reads;
* sequence homology between human and microbial genomes;
* repetitive or low-complexity sequences;
* ambiguous taxonomic assignments;
* localized cross-mapping to small genomic regions;
* reagent and laboratory contamination;
* environmental microbial background;
* PCR amplification and duplicate reads.

To address these challenges, sc-MIST applies multiple independent layers of evidence before a microbial taxon is retained.

A central design principle is that a high-confidence microbial signal should be supported by:

1. **sequence-level evidence** consistent with the assigned microbial taxon;
2. **diverse k-mer support** rather than repeated evidence from only a few sequence fragments;
3. **broad genomic support** rather than reads piling up in a highly localized region;
4. **consistent signal across cell barcodes and samples**;
5. **abundance exceeding expected contamination backgrounds**;
6. **cell-associated evidence recoverable through valid scRNA-seq barcodes and UMIs**.

---

# Workflow

```text
Raw scRNA-seq FASTQ
        |
        v
1. Human read depletion
   Cell Ranger / zUMI
        |
        v
Human-unmapped reads
        |
        v
2. Initial microbial classification
   Kraken2 / Kraken2Uniq
        |
        v
Candidate read -> taxid assignments
        |
        v
3. Sequence-level validation
   Taxon-specific BWA-MEM alignment
   - alignment length > 45 bp
   - sequence identity > 95%
        |
        v
Removal of residual human reads
   - MHC/HLA
   - human rRNA
        |
        v
4. Kraken2 re-classification
        |
        v
Clean taxonomic assignments
        |
        v
5. False-positive denoising
   - taxon-supporting k-mer fraction
   - barcode-level k-mer correlation
   - sample-level k-mer correlation
   - microbial genome read-distribution assessment
        |
        v
6. Contamination filtering
   - cell-line background filtering
   - Decontam
   - literature-derived contaminant blacklist
   - MGnify environmental blacklist
        |
        v
7. Cell-associated microbial quantification
   Cell barcode + UMI
        |
        v
High-confidence Cell x Microbe matrix
```

---

# 1. Human read depletion

sc-MIST first removes reads that can be explained by the human reference genome.

For 10x Genomics datasets, human alignment is performed using **Cell Ranger**. For compatible non-10x droplet-based datasets, barcode-aware preprocessing tools such as **zUMI** can be used.

Only reads remaining unmapped to the human reference are passed to downstream microbial analysis.

For 10x scRNA-seq data:

```text
R1 = cell barcode + UMI
R2 = cDNA sequence
```

Cell barcode and UMI information are preserved so that microbial reads can later be reassigned to individual cells.

### Human reference

In our implementation, the complete **T2T-CHM13v2.0** human reference genome was included to maximize removal of human-derived sequences.

Reference:

https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

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

---

# 2. Curated microbial reference and initial classification

## Human-associated microbial reference

To reduce false-positive assignments from environmentally abundant or poorly supported taxa, sc-MIST restricts microbial detection to taxa with prior evidence of human association.

The curated microbial reference incorporates human-associated microorganisms reported in resources including:

* MGnify
* mBodyMap
* Human Microbiome Project
* published literature

Because several cancers have established viral etiologies, the reference also includes **NCBI-annotated human-infecting viruses**.

This strategy reduces the search space while preserving microbial taxa with plausible relevance to human biology.

---

# 3. Kraken2 database construction

The sc-MIST Kraken2 database contains:

```text
Human
+
Human-associated bacteria
+
Human-infecting viruses
```

## Human reference

The human component is based on:

```text
T2T-CHM13v2.0
```

## Viral reference

Reference viral genomes were obtained from **NCBI Virus / RefSeq**.

Reference used in our analysis:

https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq

Database snapshot used in our study:

```text
2024-04-24
```

## Bacterial reference

Reference bacterial genomes can be downloaded using Kraken2:

```bash
kraken2-build --download-library bacteria --db /path/to/kraken_db
```

The downloaded bacterial collection can then be restricted to the curated human-associated taxa used by sc-MIST.

## NCBI taxonomy

Kraken2 database construction also requires the NCBI taxonomy database:

https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

Taxonomy snapshot used in our study:

```text
2024-05-27
```

## Build the Kraken2 database

After preparing the human, bacterial, viral, and taxonomy components:

```bash
kraken2-build --build --db /path/to/kraken_db
```

---

# 4. Initial microbial classification

Human-unmapped reads are classified using **Kraken2 / Kraken2Uniq**.

In our workflow, sequencing reads are decomposed into k-mers, and taxonomic assignments are inferred from k-mer matches against the custom reference database.

For Kraken2Uniq-based analyses:

```text
k = 35
```

Each read receives a putative taxonomic assignment using the lowest common ancestor framework.

Conceptually:

```text
Read
  |
  +-- k-mer 1 -> taxid A
  +-- k-mer 2 -> taxid A
  +-- k-mer 3 -> taxid B
  +-- ...
  |
  v
LCA-based read assignment
```

### Input

```text
human-unmapped FASTQ
custom Kraken2 database
```

### Output

Typical outputs include:

```text
sample.kraken.output.txt
sample.kraken.report.txt
sample.kraken.report.mpa.txt
```

The read-level output records:

```text
read ID -> assigned taxid
```

while the Kraken report summarizes microbial taxa detected at the sample level.

---

# 5. Sequence-level validation with BWA-MEM

Kraken-based taxonomic assignments are independently cross-checked using **BWA-MEM**.

Reads assigned to microbial species or genera are aligned against the corresponding taxon-specific microbial reference genome or reference library.

A read is retained only when both criteria are satisfied:

```text
alignment length > 45 bp
sequence identity > 95%
```

This additional alignment layer helps remove taxonomic assignments driven by short or weak sequence similarities.

### Example

| Read    | Kraken taxon | Alignment length | Identity | Retained |
| ------- | ------------ | ---------------: | -------: | -------- |
| Read001 | H. pylori    |            87 bp |    99.1% | Yes      |
| Read002 | H. pylori    |            36 bp |    99.0% | No       |
| Read003 | Pseudomonas  |            82 bp |    92.4% | No       |

---

# 6. Removal of residual human MHC/HLA and rRNA reads

Standard host alignment may fail to eliminate all human-derived reads, particularly sequences from highly polymorphic or repetitive regions.

Candidate microbial reads passing taxon-specific BWA validation are therefore aligned again against dedicated human references containing:

* MHC/HLA sequences;
* human ribosomal RNA;
* mitochondrial rRNA where applicable.

Reads aligning to these references are removed.

### Input

```text
BWA-validated microbial candidate reads
HLA/MHC reference
human rRNA reference
```

### Output

```text
sample.clean.R1.fastq.gz
sample.clean.R2.fastq.gz
```

These reads represent the sequence-level validated microbial candidate pool.

---

# 7. Kraken2 re-classification

Because BWA validation and host-residual filtering alter the read population, the retained reads are re-classified using Kraken2.

This step regenerates the read-level and sample-level taxonomic information required for downstream statistical denoising.

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

All downstream k-mer, read-count, and microbial-abundance calculations are based on this cleaned read population.

---

# 8. Read-level k-mer support

For each classified read, sc-MIST evaluates the fraction of informative k-mers that directly support the assigned microbial taxon.

Reads are retained when:

```text
taxon-supporting k-mer fraction >= 0.50
```

For example:

```text
18 taxon-supporting k-mers
-------------------------- = 0.90 -> retain
20 informative k-mers
```

whereas:

```text
7
-- = 0.35 -> remove
20
```

This criterion reduces assignments in which the final read-level taxonomic label is supported by only a minority of the underlying sequence evidence.

---

# 9. Barcode-level k-mer correlation

For droplet-based scRNA-seq datasets, sc-MIST evaluates microbial sequence support across cell barcodes.

For each taxon within each sample, the following quantities are calculated for each barcode:

```text
total k-mers
unique k-mers
```

A Spearman correlation is then calculated between total and unique k-mer counts across barcodes.

Taxa are retained when:

```text
Spearman r > 0.1
```

### Rationale

A genuine microbial signal should generally show increasing sequence diversity as the amount of microbial sequence evidence increases.

For example:

```text
Barcode   Total k-mers   Unique k-mers
BC01      10             3
BC02      30             8
BC03      60             17
BC04      120            31
```

is consistent with a diverse microbial signal.

In contrast:

```text
Barcode   Total k-mers   Unique k-mers
BC01      20             2
BC02      80             3
BC03      150            2
BC04      300            3
```

suggests that increasing apparent abundance is driven primarily by repeated observations of a small number of sequence fragments.

---

# 10. Sample-level k-mer correlation

Within each cohort, sc-MIST evaluates the consistency of microbial evidence across samples.

This analysis is performed when:

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
AND
r2 > 0.1
AND
r3 > 0.1
```

These three quantities capture complementary properties of microbial evidence:

```text
Reads
= number of sequencing reads assigned to the taxon

Total minimizers / k-mers
= total sequence-level taxonomic evidence

Unique minimizers / unique k-mers
= diversity of sequence evidence
```

A high-confidence microbial taxon is expected to show concordant variation among all three measurements across samples.

---

# 11. Microbial genome read-distribution assessment

sc-MIST additionally evaluates how microbial reads are distributed across the assigned microbial genome.

A plausible microbial signal should generally be supported by reads distributed across multiple genomic regions.

By contrast, reads piling up within a very narrow genomic interval can indicate:

* conserved sequence homology;
* host-microbial cross-mapping;
* repetitive sequence artifacts;
* reference-specific alignment artifacts.

Conceptually:

```text
High-confidence pattern

Genome
|------------------------------------------------|
     ^       ^           ^      ^           ^
     reads distributed across multiple regions
```

versus:

```text
Potentially spurious pattern

Genome
|------------------------------------------------|
                    ^^^^^^^
                    localized pile-up
```

Localized read distributions can therefore be flagged or removed as potential false-positive assignments.

---

# 12. Contamination filtering

Sequence-valid microbial assignments may still represent laboratory, reagent, or environmental contamination.

sc-MIST therefore applies additional contamination-control procedures.

---

## 12.1 Cell-line background filtering

For each microbial taxon, abundance is normalized as RPMM:

```text
RPMM =
taxon reads
--------------------- x 1,000,000
total microbial reads
```

The abundance observed in each biological sample is compared with the corresponding taxon-specific abundance distribution in cell-line reference datasets.

In the current workflow, the 90th percentile can be used as the background threshold:

```text
sample RPMM > cell-line q0.90 RPMM
```

Taxa failing to exceed the cell-line background are removed.

---

## 12.2 Decontam

The R package **decontam** is used as an additional contamination-detection layer.

The frequency-based method identifies potential contaminants based on their relationship with sample nucleic-acid concentration.

Contaminants often show increased relative abundance in samples with lower total nucleic-acid concentration.

The resulting contamination status can be incorporated into the final microbial filtering procedure.

---

## 12.3 Literature-derived contaminant blacklist

Taxa previously reported as common laboratory or reagent contaminants are additionally flagged using a curated literature-derived blacklist.

This provides an independent safeguard against known low-biomass sequencing contaminants.

---

## 12.4 MGnify environmental blacklist

sc-MIST also incorporates an environmental-background filter derived from MGnify.

Microorganisms repeatedly reported across multiple environment-associated studies can be considered high-risk environmental background taxa.

In our implementation, taxa reported in:

```text
> 4 independent environment-associated studies
```

can be flagged for exclusion.

---

# 13. Cell-associated microbial quantification

After false-positive and contamination filtering, retained microbial reads are linked back to individual cells using their original scRNA-seq cell barcodes.

For UMI-based droplet scRNA-seq, reads are grouped by:

```text
cell barcode
+
microbial taxon
+
UMI
```

Reads sharing the same barcode, microbial taxon, and UMI are counted once to reduce PCR amplification bias.

### Example

```text
Read001   BC01   UMI_A   H. pylori
Read002   BC01   UMI_A   H. pylori
Read003   BC01   UMI_B   H. pylori
```

Although three sequencing reads are present, the microbial abundance for this cell is:

```text
BC01 H. pylori = 2 unique UMIs
```

---

# 14. Final output

The primary sc-MIST output is a:

```text
Cell x Microbe matrix
```

Example:

| Cell barcode | H_pylori | EBV | HPV16 |
| ------------ | -------: | --: | ----: |
| BC01         |        2 |   0 |     0 |
| BC02         |        0 |   1 |     0 |
| BC03         |        4 |   0 |     2 |
| BC04         |        0 |   0 |     0 |

This matrix can be integrated directly with the corresponding host scRNA-seq expression matrix.

For example:

```text
Cell barcode
     |
     +---- host gene-expression profile
     |
     +---- cell type annotation
     |
     +---- microbial UMI counts
```

This enables analyses of:

* microbial prevalence across cell types;
* tumor-specific and cell-type-specific microbial patterns;
* microbial-positive versus microbial-negative host transcriptional states;
* inter-tumor microbial heterogeneity;
* intra-tumor microbial heterogeneity;
* host transcriptional responses associated with cell-associated microbial signals.

---

# 15. Interpretation

A microbial signal passing sc-MIST should be interpreted as:

> a high-confidence cell-associated microbial sequencing signal supported by taxonomically consistent sequence evidence and passing the implemented false-positive and contamination filters.

Detection of microbial RNA within a cell barcode does **not**, by itself, prove intracellular infection.

Possible biological origins include:

* intracellular microorganisms;
* cell-surface-associated microorganisms;
* phagocytosed microbial material;
* extracellular microbial RNA associated with a captured cell;
* droplet co-encapsulation.

Biological interpretation should therefore consider tumor context, cell identity, microbial biology, and orthogonal validation where available.

---

# 16. Input data

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

For other droplet-based platforms, the appropriate cell-barcode and UMI structure must be specified.

---

# 17. Recommended repository structure

```text
scMIST/
├── README.md
├── LICENSE
├── .gitignore
├── Snakefile
│
├── config/
│   └── config.example.yaml
│
├── workflow/
│   ├── rules/
│   └── scripts/
│
├── envs/
│   └── environment.yaml
│
├── resources/
│   ├── lab_contaminant_blacklist.tsv
│   ├── mgnify_environment_blacklist.tsv
│   └── human_associated_taxa.tsv
│
├── example/
│   ├── config.yaml
│   └── toy_data/
│
└── docs/
```

---

# 18. Configuration

Machine-specific paths should be specified through a YAML configuration file rather than hard-coded in the workflow.

Example:

```yaml
samples:
  - S01
  - S02
  - S03

fastq_dir: "/path/to/fastq"

output_dir: "results"

human_reference: "/path/to/T2T_CHM13_reference"

kraken_database: "/path/to/scMIST_kraken_database"

bwa_library: "/path/to/taxon_specific_bwa_library"

hla_reference: "/path/to/hla_reference.fa"

human_rrna_reference: "/path/to/human_rrna.fa"

cell_barcode_length: 16

umi_length: 12

bwa_min_alignment_length: 45

bwa_min_identity: 0.95

min_taxon_kmer_fraction: 0.50

barcode_correlation_threshold: 0.10

sample_correlation_threshold: 0.10

cell_line_quantile: 0.90
```

---

# 19. Software requirements

Core software used by sc-MIST includes:

```text
Snakemake
Cell Ranger
zUMI
Kraken2
Kraken2Uniq
BWA-MEM
samtools
seqkit
R
```

Relevant R packages include:

```text
data.table
tidyverse
dplyr
optparse
ShortRead
decontam
```

Software versions used for the published analysis should be recorded in the environment configuration.

---

# 20. Installation

Clone the repository:

```bash
git clone https://github.com/<USERNAME>/scMIST.git
cd scMIST
```

Create the Conda environment:

```bash
conda env create -f envs/environment.yaml
conda activate scmist
```

Software requiring separate installation, including Cell Ranger, should be installed according to the corresponding provider instructions.

---

# 21. Running sc-MIST

Example:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config/config.yaml \
    --cores 16
```

For cluster or HPC execution, an appropriate Snakemake execution profile can be supplied.

---

# 22. Example output structure

```text
results/
├── 01_human_unmapped/
├── 02_kraken_initial/
├── 03_bwa_validation/
├── 04_host_residual_filter/
├── 05_kraken_reclassification/
├── 06_false_positive_filter/
├── 07_contamination_filter/
└── 08_microbe_quantification/
```

The final microbial abundance output can be stored as:

```text
sample.microbial_UMI_matrix.tsv
```

Additional taxon-level quality-control tables are strongly recommended so that the filtering history of every detected microorganism remains traceable.

---

# 23. Methodological summary

In brief, sc-MIST identifies microbial signals using a sequence of increasingly stringent filters:

```text
Human-unmapped
      |
      v
Kraken-positive
      |
      v
Taxon-specific BWA-supported
      |
      v
MHC/rRNA-negative
      |
      v
Kraken re-confirmed
      |
      v
>=50% taxon-supporting k-mers
      |
      v
Barcode correlation r > 0.1
      |
      v
Sample correlations
r1 > 0.1
r2 > 0.1
r3 > 0.1
      |
      v
Genome-distribution support
      |
      v
Cell-line / Decontam / blacklist filtering
      |
      v
High-confidence cell-associated microorganism
```

---

# 24. Citation

sc-MIST was developed for:

**A Pan-Cancer Cell-Associated Microbiome Atlas Reveals Inter- and Intra-Tumor Heterogeneity**
*Under Review*

A formal citation will be added upon publication.

If you use sc-MIST before publication, please cite this GitHub repository and acknowledge the manuscript above.

---

# 25. Acknowledgements

sc-MIST incorporates concepts from prior work on k-mer-based microbial denoising in single-cell sequencing, including the SAHMI framework. We thank the developers of SAHMI and the broader open-source microbial genomics community for providing useful methodological foundations and software resources.

---

# 26. License

This project is distributed under the **MIT License**.

See:

```text
LICENSE
```

for details.

---

# 27. Contact

For questions, bug reports, or feature requests, please open an issue in this GitHub repository.
