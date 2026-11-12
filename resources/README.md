# Bundled resources

These files are small enough to distribute. They are used at the steps shown
below; large reference databases are **not** included.

| File | What it contains | Used in |
|---|---|---|
| `curated_bacteria_species.tsv` | 3,808 bacterial species retained after scoring human-association evidence from PMID 32467386, HMP, MGnify, and mBodyMap. Columns include species/genus/phylum names and taxids, score, and evidence at each rank. | Step 04 |
| `kraken_overlap_databases_score6.xlsx` | The 3,808 species with score >= 6, reduced to seven columns: `species`, `species_taxid`, `genus`, `genus_taxid`, `phylum`, `phylum_taxid`, `score`. | Reference audit |
| `CCLE_quantile_mhc_rrna.txt` | Cell-line RPMM quantile table. Step 08 reads `taxid` and `threshold_0.9_exclude0`; HBV taxid `10407` is overridden to `10.507194`. | Step 08 |
| `HVIDB_PPIs.txt` | HVIDB human-virus protein interaction table. Step 11 matches `Organism_Interactor_virus` (viral taxid). | Step 11 |
| `humanvirusbase_Virus_Organ_pairs.csv` | Human Virus Database organ-pair table. Step 11 matches `VirusName`. | Step 11 |
| `black_list.txt` | Laboratory/reagent contaminant genera from the published contamination list (Salter et al., *BMC Biology* 2014, Table 1: https://doi.org/10.1186/s12915-014-0087-z). | Step 10 |
| `mgnify_black.txt` | MGnify genera reported in more than four environment-associated studies (`count_n > 4`). Source: MGnify, https://www.ebi.ac.uk/metagenomics/. | Step 10 |
| `ccle_blacklist_taxids.txt` | Taxids excluded after the cell-line RPMM screen against the public cell-line microbiome data used by De et al. | Step 10 |

## Scoring used for `curated_bacteria_species.tsv`

For every bacterial species downloaded for the Kraken2 library, a hit in each
external resource was scored by the taxonomic level at which the hit occurred:

```text
species = 3
genus   = 2
phylum  = 1
```

Scores were summed across PMID 32467386, HMP, MGnify, and mBodyMap. Species
with a total score **>= 6** were retained. The final set contains 3,808
species. The compact workbook contains only the rows and columns needed to
audit this threshold; it does not contain all original evidence columns.

## Viral host evidence and the NCBI metadata file

Step 11 keeps a virus only when at least two of NCBI Virus, HVD, and HVIDB
support a human association. HVD and HVIDB are bundled here.

The NCBI Virus nucleotide host metadata export used for the third evidence
source is about **1.3 GB** and is therefore not committed. Download it from:

```text
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq
```

The file must contain these columns:

```text
Organism_Name  Species  Genus  Host
```

Then set `VIRUS_NCBI_HOST` in `config/config.sh`. If this file is unavailable,
the script can still use HVD and HVIDB, but the evidence score will be based on
fewer sources; document that change in your analysis.

## Not bundled

The Kraken2 database, taxon-specific BWA indices, HLA/rRNA FASTA, and taxonkit
lineage matrix are too large or are derived resources. See
[`docs/01_reference_build.md`](../docs/01_reference_build.md) and the main
README for construction instructions.
