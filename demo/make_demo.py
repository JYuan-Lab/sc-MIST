#!/usr/bin/env python3
"""Create a tiny, fully synthetic 3-sample sc-MIST demo dataset.

The generated data are intentionally small enough to commit to GitHub. They
exercise every step of the pipeline and are not intended for biological
interpretation.
"""

from __future__ import annotations

import gzip
import random
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "demo"
REF = DEMO / "reference"
INPUT = DEMO / "input"
RNG = random.Random(20241111)

TAXA = {
    "1001": ("ToyspeciesA", "1000", "ToygenusA", "100", "ToyphylumA", "Bacteria"),
    "1002": ("ToyspeciesA2", "1000", "ToygenusA", "100", "ToyphylumA", "Bacteria"),
    "2001": ("ToyspeciesB", "2000", "ToygenusB", "200", "ToyphylumB", "Bacteria"),
    "3001": ("Toyvirus", "3000", "Toyvirusgenus", "", "", "Viruses"),
    "9606": ("Homo_sapiens", "", "", "", "", "Eukaryota"),
}


def dna(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main() -> None:
    for d in (REF, INPUT):
        d.mkdir(parents=True, exist_ok=True)

    shared = dna(250)
    refs = {
        "1001": shared + dna(250),
        "1002": shared + dna(250),
        "2001": dna(500),
        "3001": dna(500),
        "9606": dna(500),
    }

    # Taxon-specific references used by Kraken2. Genus-level BWA references
    # are built later by concatenating their species members.
    for taxid, seq in refs.items():
        write(REF / f"taxon_{taxid}.fa", f">kraken:taxid|{taxid}\n{seq}\n")

    # Host residual references: one sequence resembles a 1001 fragment, so a
    # small number of microbial-looking reads can be removed in step 05.
    hla = refs["1001"][250:450] + dna(100)
    rrna = dna(300)
    write(REF / "host_residual.fa", f">toy_HLA\n{hla}\n>toy_rRNA\n{rrna}\n")

    taxonomy = REF / "kraken_taxonomy"
    taxonomy.mkdir(parents=True, exist_ok=True)
    nodes = [
        ("1", "1", "no rank"),
        ("2", "1", "superkingdom"),
        ("10239", "1", "superkingdom"),
        ("100", "2", "phylum"),
        ("200", "2", "phylum"),
        ("1000", "100", "genus"),
        ("2000", "200", "genus"),
        ("3000", "10239", "genus"),
        ("1001", "1000", "species"),
        ("1002", "1000", "species"),
        ("2001", "2000", "species"),
        ("3001", "3000", "species"),
        ("9606", "1", "species"),
    ]
    nodes_text = "".join(
        f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
        for t, p, r in nodes
    )
    write(taxonomy / "nodes.dmp", nodes_text)

    names = [
        ("1", "root"),
        ("2", "Bacteria"),
        ("10239", "Viruses"),
        ("100", "ToyphylumA"),
        ("200", "ToyphylumB"),
        ("1000", "ToygenusA"),
        ("2000", "ToygenusB"),
        ("3000", "Toyvirusgenus"),
        ("1001", "ToyspeciesA"),
        ("1002", "ToyspeciesA2"),
        ("2001", "ToyspeciesB"),
        ("3001", "Toyvirus"),
        ("9606", "Homo sapiens"),
    ]
    names_text = "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in names)
    write(taxonomy / "names.dmp", names_text)

    # Lineage matrix consumed by step 04.
    lineage = [
        "taxid\tkindom\tphylum\tgenus\tspecies\tg_taxid\ts_taxid",
        "1000\tk__Bacteria\tp__ToyphylumA\tg__ToygenusA\tNA\t1000\tNA",
        "1001\tk__Bacteria\tp__ToyphylumA\tg__ToygenusA\ts__ToyspeciesA\t1000\t1001",
        "1002\tk__Bacteria\tp__ToyphylumA\tg__ToygenusA\ts__ToyspeciesA2\t1000\t1002",
        "2000\tk__Bacteria\tp__ToyphylumB\tg__ToygenusB\tNA\t2000\tNA",
        "2001\tk__Bacteria\tp__ToyphylumB\tg__ToygenusB\ts__ToyspeciesB\t2000\t2001",
        "3000\tk__Viruses\tNA\tg__Toyvirusgenus\tNA\t3000\tNA",
        "3001\tk__Viruses\tNA\tg__Toyvirusgenus\ts__Toyvirus\t3000\t3001",
        "9606\tk__Eukaryota\tNA\tNA\ts__Homo_sapiens\tNA\t9606",
    ]
    write(REF / "taxonkit_lineage_matrix.tsv", "\n".join(lineage) + "\n")

    curated = [
        "species\tspecies_taxid\tgenus\tgenus_taxid\tphylum\tphylum_taxid\tscore",
        "s__ToyspeciesA\t1001\tg__ToygenusA\t1000\tp__ToyphylumA\t100\t6",
        "s__ToyspeciesA2\t1002\tg__ToygenusA\t1000\tp__ToyphylumA\t100\t6",
        "s__ToyspeciesB\t2001\tg__ToygenusB\t2000\tp__ToyphylumB\t200\t6",
    ]
    write(REF / "curated_bacteria_species.tsv", "\n".join(curated) + "\n")

    q90 = ["taxid\tq90", "1000\t0", "1001\t0", "1002\t0", "2001\t0", "3001\t0"]
    write(REF / "cell_line_rpmm_q90.tsv", "\n".join(q90) + "\n")
    write(REF / "black_list_demo.txt", "phylum\tgenus\n")
    write(
        REF / "mgnify_black_demo.csv",
        "taxonomy_degree,taxonomy,env1,env2,count_n\n",
    )
    write(REF / "ccle_blacklist_demo.txt", "999999\n")

    barcodes = [
        "AAACCCGGGTTTAAAA",
        "CCCGGGTTTAAAACCC",
        "GGGTTTAAAACCCGGG",
        "TTTAAAACCCGGGTTT",
        "ACGTACGTACGTACGT",
        "TGCATGCATGCATGCA",
        "GATTACAGATTACAGA",
        "CTAGCTAGCTAGCTAG",
    ]

    samples = ["sample1", "sample2", "sample3"]
    # Per-sample read multiplier; each barcode gets 1-4 copies of each taxon.
    multipliers = [1, 2, 3]
    taxon_order = ["1000", "1001", "1002", "2001", "3001"]

    for si, sample in enumerate(samples):
        sample_dir = INPUT / sample / "outs"
        matrix_dir = sample_dir / "filtered_feature_bc_matrix"
        matrix_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "possorted_genome_bam.bam").parent.mkdir(parents=True, exist_ok=True)

        with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt") as fh:
            for bc in barcodes:
                fh.write(bc + "-1\n")
        with gzip.open(matrix_dir / "features.tsv.gz", "wt") as fh:
            for i in range(1, 11):
                fh.write(f"toyGene{i}\ttoyGene{i}\tGene Expression\n")
        with gzip.open(matrix_dir / "matrix.mtx.gz", "wt") as fh:
            fh.write("%%MatrixMarket matrix coordinate integer general\n")
            fh.write("%\n")
            fh.write(f"10 8 80\n")
            for g in range(1, 11):
                for b in range(1, 9):
                    fh.write(f"{g} {b} 1\n")

        sam = sample_dir / "possorted_genome_bam.sam"
        with sam.open("w") as fh:
            fh.write("@HD\tVN:1.6\tSO:unsorted\n")
            fh.write("@SQ\tSN:chr1\tLN:1000000\n")
            read_no = 0
            for ti, taxid in enumerate(taxon_order):
                for bi, bc in enumerate(barcodes):
                    n_reads = multipliers[si] * (1 + ((bi + ti) % 4))
                    for copy in range(n_reads):
                        read_no += 1
                        umi = dna(10)
                        rid = f"{sample}_read{read_no:05d}"
                        if taxid == "1000":
                            start = (copy * 17 + bi * 5) % 100
                            seq2 = refs["1001"][start:start + 150]
                        else:
                            start = 250 + ((copy * 17 + bi * 5) % 100)
                            seq2 = refs[taxid][start:start + 150]
                        seq1 = bc + umi + dna(124)
                        fh.write(
                            f"{rid}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq1}\t{'I' * len(seq1)}\n"
                        )
                        fh.write(
                            f"{rid}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq2}\t{'I' * len(seq2)}\n"
                        )

            # A few human-classified reads and host-residual reads demonstrate
            # that the host filters are active.
            for _ in range(2):
                read_no += 1
                rid = f"{sample}_host{read_no:05d}"
                seq1 = barcodes[0] + dna(10) + dna(124)
                seq2 = refs["9606"][RNG.randrange(0, 350):][:150]
                fh.write(f"{rid}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq1}\t{'I' * len(seq1)}\n")
                fh.write(f"{rid}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq2}\t{'I' * len(seq2)}\n")
            for _ in range(1):
                read_no += 1
                rid = f"{sample}_residual{read_no:05d}"
                seq1 = barcodes[1] + dna(10) + dna(124)
                seq2 = refs["1001"][260:410]
                fh.write(f"{rid}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq1}\t{'I' * len(seq1)}\n")
                fh.write(f"{rid}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq2}\t{'I' * len(seq2)}\n")

        subprocess.run(
            ["samtools", "view", "-b", "-o", str(sample_dir / "possorted_genome_bam.bam"), str(sam)],
            check=True,
        )
        sam.unlink()

    write(
        DEMO / "samples.tsv",
        "sample_id\tcellranger_outs\n"
        "sample1\tdemo/input/sample1/outs\n"
        "sample2\tdemo/input/sample2/outs\n"
        "sample3\tdemo/input/sample3/outs\n",
    )
    print("Demo input and source references written under", DEMO)


if __name__ == "__main__":
    main()
