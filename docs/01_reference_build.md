# Reference construction

The pipeline needs several external resources. Only the small, auditable tables
are in `resources/`; the large sequence databases and NCBI Virus metadata must
be downloaded or built once and then reused across projects.

```text
reference/
├── kraken_db/                      # KRAKEN_DB
├── bwa_library/                    # BWA_LIBRARY
├── host_residual.fa                # HOST_RESIDUAL_FASTA (+ BWA index)
├── taxonkit_lineage_matrix.tsv     # TAXON_LINEAGE
└── cell_line_rpmm_q90.tsv          # CELLLINE_Q90
```

> **Time and disk.** Building the full Kraken2 database and the taxon-specific
> BWA indices is the expensive part of this project. Plan for tens to hundreds
> of GB of free disk, depending on how many bacterial assemblies you retain.
> The commands below are written for one-time reference preparation, not for
> every project.

After building the references, check them with:

```bash
cp config/config.example.sh config/config.sh
$EDITOR config/config.sh
bash scripts/utils/check_references.sh --config config/config.sh
```

The checker only validates paths, headers, and the presence of indices; it does
not validate biological content.

## 1. Kraken2 database

### 1.1 Contents

The database used in the manuscript contains four components:

```text
Human: T2T-CHM13v2.0
Bacteria: Kraken2 bacterial reference genomes, restricted to taxa with prior evidence of human association
Viruses: human-infecting viral genomes from NCBI Virus / RefSeq
Taxonomy: NCBI taxonomy
```

The source snapshots used for the released analysis were:

| Component | Source | Snapshot / note |
|---|---|---|
| Human genome | T2T-CHM13v2.0, `chm13v2.0.fa.gz` | `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz` |
| Bacteria | Kraken2 `bacteria` library, filtered with `resources/curated_bacteria_species.tsv` | 3,808 species; score >= 6 |
| Viruses | NCBI Virus / RefSeq, human host filter | download date 2024-04-24 |
| Taxonomy | NCBI `taxdump.tar.gz` | 2024-05-27 |
| Kraken2 | bioconda Kraken2 | v2.1.3; k-mer length 35 |

The bacterial species universe was defined by comparing Kraken2 bacterial taxa
with PMID 32467386, HMP, MGnify, and mBodyMap. A hit at species, genus, or
phylum was scored 3, 2, or 1, respectively; species with a total score >= 6
were retained. The final list is:

```text
resources/curated_bacteria_species.tsv
```

The compact audit table is:

```text
resources/kraken_overlap_databases_score6.xlsx
```

It contains only the seven essential columns and the 3,808 species with
`score >= 6`. If you want a stricter `score > 6` set, filter the `score`
column again.

### 1.2 Create the database and add the taxonomy

```bash
export DB=/path/to/reference/kraken_db
export STAGING=/path/to/reference/kraken_staging
mkdir -p "$DB" "$STAGING"

# Option A: let Kraken2 download the current taxonomy.
kraken2-build --download-taxonomy --db "$DB"

# Option B: pin the manuscript snapshot (recommended for exact reproduction).
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
mkdir -p "$DB/taxonomy"
cp nodes.dmp names.dmp "$DB/taxonomy/"
```

If you use Option B, also copy `merged.dmp` and `delnodes.dmp` if they are
present; taxonkit can use them to resolve merged or deleted taxids.

### 1.3 Add human T2T-CHM13v2.0

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
kraken2-build --add-to-library chm13v2.0.fa.gz --db "$DB"
```

Kraken2 accepts gzipped FASTA directly. The human sequence must be present in
the database because the first classification pass is intended to classify
host reads as human rather than as microbial.

### 1.4 Add the curated bacterial genomes

Download the Kraken2 bacterial library into a staging database:

```bash
kraken2-build --download-library bacteria --db "$STAGING"
```

This downloads many assemblies and can take a long time. The downloaded
`assembly_summary.txt` has one row per assembly. Column 1 is the accession,
column 6 is the taxid, and column 7 is the species taxid. Keep assemblies whose
species taxid is in the curated list:

```bash
tail -n +2 resources/curated_bacteria_species.tsv \
  | cut -f2 \
  | sort -u > species_taxids.txt

awk -F'\t' 'NR == FNR {keep[$1]; next}
            !/^#/ && (($6 in keep) || ($7 in keep)) {print $1}' \
  species_taxids.txt "$STAGING/library/bacteria/assembly_summary.txt" \
  > selected_assembly_accessions.txt

wc -l selected_assembly_accessions.txt
```

Add the selected assemblies to the final database:

```bash
while read -r acc; do
  f=$(find "$STAGING/library/bacteria" -maxdepth 1 -name "${acc}*_genomic.fna.gz" | head -n 1)
  if [[ -z "$f" ]]; then
    echo "WARNING: no FASTA for $acc" >&2
  else
    kraken2-build --add-to-library "$f" --db "$DB"
  fi
done < selected_assembly_accessions.txt
```

If the Kraken2 release you use stores files under a different name, inspect
`$STAGING/library/bacteria` and adapt the `find` pattern. If several assemblies
exist for one species, either keep all of them or choose one representative
assembly per species and record the accessions in your methods. Do not add the
entire unfiltered bacterial library if you want the same species universe as
the manuscript.

### 1.5 Add human-infecting viruses

Download viral genomes from NCBI Virus / RefSeq with a human host filter:

```text
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq
```

The released analysis used the 2024-04-24 download. The exact accession list is
not bundled because NCBI Virus records and host annotations change over time.
Apply the following rule:

1. keep records annotated as human-infecting in NCBI Virus;
2. cross-check with HVD (Human Virus Database; PMID 35435220) and HVIDB
   (PMID 33515030), and require support from at least two of the three sources
   (NCBI Virus, HVD, HVIDB);
3. export the genomes as FASTA and add them to the same Kraken2 database.

```bash
kraken2-build --add-to-library human_infecting_viruses.fa --db "$DB"
```

Every viral FASTA header must either contain a Kraken2 taxid
(`>kraken:taxid|<taxid>`) or be resolvable through a `seqid2taxid.map` file.
Check a few headers before building:

```bash
grep '^>' human_infecting_viruses.fa | head
```

### 1.6 Build and check

```bash
kraken2-build --build --db "$DB" --threads 8 --kmer-len 35
ls "$DB"/{hash.k2d,opts.k2d,taxo.k2d}
```

The manuscript used Kraken2 v2.1.3, `k-mer = 35`, and the snapshots listed
above. The pipeline can read a database built by either `kraken2` or
`Kraken2Uniq`; set `KRAKEN2_BIN` in `config/config.sh` to the executable you
used. `Kraken2Uniq` was used in the released analysis.

Set:

```bash
KRAKEN_DB=/path/to/reference/kraken_db
```

## 2. Taxon-specific BWA library

Step 04 validates each Kraken2 assignment by aligning the read to a reference
for the assigned taxid. The library is a directory of BWA-indexed FASTA files
named by taxid:

```text
bwa_library/
├── 1001.fa
├── 1001.fa.amb
├── 1001.fa.ann
├── 1001.fa.bwt
├── 1001.fa.pac
├── 1001.fa.sa
├── 1000.fa               # genus = concatenation of its species
└── ...
```

Rules:

- **Species taxid:** use the species genome/sequence FASTA.
- **Genus taxid:** concatenate the FASTAs of all retained species in that
  genus, then index the combined FASTA.
- Create one entry for every taxid that step 04 may receive: bacterial species,
  bacterial genera, and viral species/genera present in the Kraken2 database.
- Run `bwa index` on every FASTA.

### 2.1 Prepare one FASTA per species taxid

Create a directory of files named `<taxid>.fa`, `<taxid>.fasta`, or
`<taxid>.fna`. The FASTA can come from the same genomes used to build the
Kraken2 database. If a downloaded file contains several contigs for one
species, that is fine: keep them in one FASTA with multiple `>` records.

### 2.2 Build the lineage matrix first

The BWA helper needs the seven-column lineage table described in section 4.
Build a taxid list from the curated bacteria table (species and genus) and add
the viral taxids you retained:

```bash
tail -n +2 resources/curated_bacteria_species.tsv \
  | awk -F'\t' '{print $2; print $4}' \
  | sort -u > reference_taxids.txt

# Add viral taxids, one per line, for example:
# cat virus_taxids.txt >> reference_taxids.txt
```

Then build the lineage matrix:

```bash
Rscript scripts/utils/build_lineage_matrix.R \
  --taxids reference_taxids.txt \
  --taxdump /path/to/taxdump \
  --out /path/to/reference/taxonkit_lineage_matrix.tsv
```

See section 4 for the required columns and an alternative using an existing
taxonomy table.

### 2.3 Build the BWA library

```bash
bash scripts/utils/build_bwa_library.sh \
  --species-dir /path/to/species_fasta \
  --lineage /path/to/reference/taxonkit_lineage_matrix.tsv \
  --outdir /path/to/reference/bwa_library
```

The helper copies each species FASTA to `<taxid>.fa`, indexes it, and creates
genus-level FASTA files by concatenating the species that share `g_taxid`. It
does not create references for taxids that have no species FASTA; check the
warnings and add the missing files if needed.

Set:

```bash
BWA_LIBRARY=/path/to/reference/bwa_library
```

## 3. HLA + human rRNA reference

Residual reads from the highly polymorphic HLA/MHC region and from human
ribosomal RNA are removed in step 05. The pipeline expects one combined,
BWA-indexed FASTA:

```bash
cat hla_gen.fasta ncbi_human_18_28S_MT_rrna.fasta > host_residual.fa
bwa index host_residual.fa
```

Sources:

- **HLA:** IPD-IMGT/HLA, `ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/`,
  file `hla_gen.fasta`.
- **Human rRNA:** an NCBI-derived FASTA containing **18S rRNA, 28S rRNA, and
  mitochondrial rRNA**. The released file was
  `ncbi_human_18_28S_MT_rrna.fasta`.

The combined reference contains:

- **HLA/MHC** sequences from IPD-IMGT/HLA;
- **Human rRNA reference**, including:
  - 18S rRNA;
  - 28S rRNA;
  - mitochondrial rRNA.

Set:

```bash
HOST_RESIDUAL_FASTA=/path/to/reference/host_residual.fa
```

Step 05 checks for `${HOST_RESIDUAL_FASTA}.bwt` and stops with a clear error if
the BWA index is missing.

## 4. Taxonkit lineage matrix

Step 04 uses a TSV with exactly these columns:

```text
taxid	kindom	phylum	genus	species	g_taxid	s_taxid
```

- `kindom` is intentionally spelled this way to match the released matrix.
- `g_taxid` is the genus taxid; `s_taxid` is the species taxid.
- Use `NA` when a rank is absent.
- Rows for both species and genus taxids are allowed.

Example:

```text
taxid	kindom	phylum	genus	species	g_taxid	s_taxid
1000	k__Bacteria	p__ToyphylumA	g__ToygenusA	NA	1000	NA
1001	k__Bacteria	p__ToyphylumA	g__ToygenusA	s__ToyspeciesA	1000	1001
1002	k__Bacteria	p__ToyphylumA	g__ToygenusA	s__ToyspeciesA2	1000	1002
3001	k__Viruses	NA	g__Toyvirusgenus	s__Toyvirus	3000	3001
```

The released matrix is ~302 MB and is not committed. Rebuild it from the same
taxdump used for the Kraken2 database:

```bash
Rscript scripts/utils/build_lineage_matrix.R \
  --taxids reference_taxids.txt \
  --taxdump /path/to/taxdump \
  --out /path/to/reference/taxonkit_lineage_matrix.tsv
```

The helper runs `taxonkit lineage -t -R`, keeps the superkingdom/kingdom,
phylum, genus, and species ranks, prefixes the names with `k__`, `p__`, `g__`,
and `s__`, and writes the exact seven columns above.

If you already have a lineage table with different column names, rename the
columns to the seven names above; that is the only hard requirement.

Set:

```bash
TAXON_LINEAGE=/path/to/reference/taxonkit_lineage_matrix.tsv
```

## 5. Cell-line RPMM 90th-percentile table

Step 08 compares each taxon's RPMM in the sample with the 90th percentile of
nonzero RPMM values observed across cell-line controls. Required columns are:

```text
taxid	q90
```

or the original column names:

```text
taxid	rank	threshold_0.9_exclude0
```

The script accepts either `q90` or `threshold_0.9_exclude0`. HBV (taxid
`10407`) is overridden to `10.507194`, as in the released analysis.

The released table is bundled as `resources/CCLE_quantile_mhc_rrna.txt`. It was
derived from the public cell-line microbiome data used by De et al. The raw
file is available from the link in the manuscript:

```text
https://www.dropbox.com/s/r6xvw1589lqyqts/cell.lines.txt?dl=0
```

The calculation is the 90th percentile of **nonzero** RPMM entries for each
taxon, where RPMM is taxon counts divided by total bacterial and viral counts
times 1,000,000. If your cell-line data are a matrix with one row per taxid and
one column per cell line, this R snippet produces the table:

```r
mat <- read.delim("cell.lines.rpmm.tsv", row.names = 1, check.names = FALSE)
q90 <- apply(mat, 1, function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) NA_real_ else as.numeric(quantile(v, 0.9, names = FALSE))
})
out <- data.frame(taxid = rownames(mat), q90 = q90)
write.table(out, "cell_line_rpmm_q90.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
```

If you have long-format data with columns `taxid`, `cell_line`, and `rpmm`,
group by `taxid` first and apply the same `quantile(..., 0.9)` to the positive
values.

If the table is unavailable, a table with `q90 = 0` for all taxa disables the
cell-line filter. That is less stringent and must be stated explicitly in your
methods.

Set:

```bash
CELLLINE_Q90="$PROJECT_ROOT/resources/CCLE_quantile_mhc_rrna.txt"
```

## 6. Blacklists

All three are bundled and used by step 10:

```text
resources/black_list.txt
resources/mgnify_black.txt
resources/ccle_blacklist_taxids.txt
```

- `black_list.txt`: laboratory/reagent contaminant genera from the published
  contamination list. It has a `phylum`/`genus` header.
- `mgnify_black.txt`: MGnify genera reported in more than four
  environment-associated studies. It has the MGnify export header.
- `ccle_blacklist_taxids.txt`: taxids excluded after the cell-line screen; one
  taxid per line, no header.

Keep the headers unchanged. If you update a list, record the source and date in
your analysis notebook.

## Viral host evidence files

Step 11 needs three evidence sources for viral retention:

```text
NCBI Virus host metadata   VIRUS_NCBI_HOST   ~1.3 GB, not bundled
HVD                        VIRUS_HVD         bundled
HVIDB                      VIRUS_HVIDB       bundled
```

The NCBI metadata download is:

```text
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq
```

The required columns are:

```text
Organism_Name  Species  Genus  Host
```

The default rule keeps a virus when at least two of the three sources support
a human association (`VIRUS_MIN_SCORE=2`).

## Final reference checklist

```text
[ ] KRAKEN_DB/hash.k2d, opts.k2d, taxo.k2d exist
[ ] BWA_LIBRARY/<taxid>.fa.bwt exists for every taxid used in step 04
[ ] HOST_RESIDUAL_FASTA and HOST_RESIDUAL_FASTA.bwt exist
[ ] TAXON_LINEAGE has the seven required columns
[ ] CELLLINE_Q90 has taxid + q90 (or threshold_0.9_exclude0)
[ ] the three blacklist files exist and keep their headers
[ ] VIRUS_NCBI_HOST, VIRUS_HVD, and VIRUS_HVIDB are readable
[ ] bash scripts/utils/check_references.sh --config config/config.sh passes
```
