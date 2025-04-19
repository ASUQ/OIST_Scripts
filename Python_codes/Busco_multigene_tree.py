#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    Busco_multigene_tree.py \
        -i <INPUT_DIR> -f <FRACTION> [-c <CORES>] \
        [-m <MAFFT_OPTION>] [-t <TRIMAL_OPTION>] [-a <AMAS_OPTION>] [-q <IQTREE_OPTION>] \
        [-o <OUT_DIR>] [-h|--help]

Purpose:
    Create multi-gene phylogenomic tree from BUSCO single-copy orthologs.

Options:
    -h, --help                  Show this
    -i, --input_dir             Path to the input directory which contains all BUSCO outputs
    -f, --fraction              Comma-spliced fractions for creating mulit-gene phylogenetic tree
    -c, --cores                 Number of CPUs to execute [default: 8]
    -m, --mafft                 mafft options [default: --globalpair --maxiterate 1000 --thread $CORES]
    -t, --trimal                trimAl options
    -a, --amas                  AMAS options
    -q, --iqtree                IQ-TREE options
    -o, --out_dir               Output directory path [default: current directory]

Required Packages and Softwares:
    Biopython, mafft, trimAl, AMAS, IQ-TREE

Author: Arno Hagenbeek (ArnoHagenbeek), Akito Shima (ASUQ)
Email: akito.shima@oist.jp
"""

import argparse
import logging
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Build multi-gene phylogenomic tree from BUSCO single-copy orthologous gene",
        epilog="Required package and Softwares: Biopython, mafft, trimAl, AMAS, IQ-TREE",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--input_dir",
        type=Path,
        required=True,
        help="Path to the input directory which contains all BUSCO outputs",
    )

    required.add_argument(
        "-f",
        "--fraction",
        type=str,
        required=True,
        help="Comma-spliced fractions for creating mulit-gene phylogenetic tree.",
    )

    optional.add_argument(
        "-c",
        "--cores",
        type=int,
        default=8,
        help="Number of CPUs to execute [default: 8]",
    )

    optional.add_argument(
        "-m",
        "--mafft",
        nargs="+",
        metavar="MAFFT_OPTION",
        default=["--globalpair", "--maxiterate", "1000", "--thread", "$CORES"],
        help="Command option for mafft [default: mafft --globalpair --maxiterate 1000 --thread $CORES]",
    )

    optional.add_argument(
        "-t",
        "--trimal",
        nargs="+",
        metavar="TRIMAL_OPTION",
        default=["--automated1"],
        help="Command option for mafft [default: trimal --automated1]",
    )

    optional.add_argument(
        "-a",
        "--amas",
        nargs="+",
        metavar="AMAS_OPTION",
        default=[""],
        help="Command option for AMAS [default: ]",
    )

    optional.add_argument(
        "-q",
        "--iqtree",
        nargs="+",
        metavar="IQTREE_OPTION",
        default=["-B", "1000", "-alrt", "1000", "-T", "$CORES"],
        help="Command option for IQ-TREE [default: ]",
    )

    optional.add_argument(
        "-o",
        "--out_dir",
        type=Path,
        default=Path("./"),
        help="Output directory path  [default: ./]",
    )

    return parser.parse_args()


def run_cmd(cmd: list[str], stdout=None) -> None:
    """Execute a subprocess call, abort on failure"""
    try:
        subprocess.check_call(cmd, stdout=stdout)
    except subprocess.CalledProcessError as e:
        logging.fatal(f"Error: command failed: {' '.join(cmd)}")
        sys.exit(e.returncode)


def collect_gene_seqs(
    input_dir: Path, output_dir: Path
) -> tuple[dict[str, set[str]], set[str]]:
    """
    Parse single-copy BUSCO FASTAs with SeqIO, prefix record.id with organism name,
    and append into per-gene files.
    """

    # gene_name: {set of orgs that have it}
    gene_dict: dict[str, set[str]] = defaultdict(set)
    org_set: set[str] = set()  # all organism names seen

    for seq_dir in input_dir.rglob("single_copy_busco_sequences"):
        if not seq_dir.is_dir():
            continue

        ###
        # Extract sample name
        # assuming structure: /input_dir/subdirectories/sample/busco_output/run_lineage/single_copy_busco_sequences/*.faa
        ###
        org_name = seq_dir.relative_to(input_dir).parts[-4]
        org_set.add(org_name)

        for faa_file in seq_dir.glob("*.faa"):
            gene = faa_file.stem
            gene_dict[gene].add(org_name)
            out_file = output_dir / f"{gene}.faa"
            with faa_file.open("r") as inp, out_file.open("a") as out:
                for rec in SeqIO.parse(inp, "fasta"):
                    rec.id = org_name
                    rec.description = ""
                    SeqIO.write(rec, out_file, "fasta")

    return gene_dict, org_set


def select_shared_genes(
    gene_dict: dict[str, set[str]], org_set: set[str], fractions: list[float]
) -> dict[float, list[str]]:
    """For each fraction, return list of genes present above the threshold.

    Args:
        gene_dict (dict): gene -> set(orgs)
        org_set (set): set of all organisms
        fractions (list): array of fractions

    Returns:
        dict: fraction -> tuple of genes present above the threshold
    """

    total = len(org_set)
    fract_dict: dict[float, list[str]] = {}
    for frac in fractions:
        threshold = total * frac
        fract_dict[frac] = [
            gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold
        ]
    return fract_dict


def write_gene_lists(fract_dict: dict, output_dir: Path) -> None:
    """Write out which genes are analyzed per fraction."""

    for frac, genes in fract_dict.items():
        file_path = output_dir / f"frac{int(frac * 100)}%_analyzed_genes.txt"
        with file_path.open("w") as out:
            out.write(f"Number of genes considered: {len(genes)}\n")
            out.write("Analyzed genes:\n")
            out.write("\n".join(genes) + "\n")


def align_and_trim(
    fract_dict: dict[float, list[str]],
    output_dir: Path,
    mafft_opts: list[str],
    trimal_opts: list[str],
) -> None:
    """Align and trim all genes in the most inclusive set

    Args:
        fract_dict (dict): fraction -> tuple of genes present above the threshold
        output_dir (Path): Path of output directory
        mafft_opts (tuple): mafft command set
        trimal_opts (tuple): trimal command set
    """
    smallest_frac = min(fract_dict.keys())
    genes = fract_dict[smallest_frac]

    for gene in genes:
        infile = output_dir / f"{gene}.faa"
        aligned = output_dir / f"{gene}_aligned.faa"
        trimmed = output_dir / f"{gene}_trimmed.faa"

        logging.info(f"Running mafft on {infile}")
        with aligned.open("w") as out_f:
            run_cmd(["mafft"] + mafft_opts + [str(infile)], stdout=out_f)

        logging.info(f"Running trimAl on {aligned}")
        run_cmd(["trimal"] + trimal_opts + ["-in", str(aligned), "-out", str(trimmed)])


def concat_alignments(
    fract_dict: dict[float, list[str]], output_dir: Path
) -> dict[float, Path]:
    """For each fraction, concatenate trimmed gene alignments"""

    cafiles: dict[float, Path] = {}
    for frac, genes in fract_dict.items():
        # concat: organism -> gene -> sequence
        concat: dict[str, dict[str, str]] = defaultdict(lambda: defaultdict(str))
        for gene in genes:
            trimmed_faa_path = output_dir / f"{gene}_trimmed.faa"
            for rec in SeqIO.parse(trimmed_faa_path, "fasta"):
                org = rec.id
                concat[org][gene] = str(rec.seq)

        # pad missing genes
        ref = next(iter(concat.values()))
        for org, seqs in concat.items():
            for gene in genes:
                if gene not in seqs:
                    seqs[gene] = "-" * len(ref[gene])

        # Write concatenated file
        cafile = output_dir / f"frac{int(frac * 100)}%_concat.faa"
        records: list[SeqRecord] = []
        for org, seqs in concat.items():
            full = "".join(seqs[gene] for gene in genes)
            records.append(SeqRecord(Seq(full), id=org, description=""))
        SeqIO.write(records, cafile, "fasta")
        cafiles[frac] = cafile
    return cafiles


def run_iqtree(cafiles: dict[float, Path], iqtree_opts: list[str]) -> None:
    """Run IQ-TREE"""

    for frac, cafile in cafiles.items():
        results_dir = cafile.parent / f"frac{int(frac * 100)}%_results"
        results_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"Running IQ-TREE on {cafile}")
        run_cmd(
            ["iqtree"]
            + iqtree_opts
            + [
                "-s",
                str(cafile),
                "-pre",
                str(results_dir / cafile.name),
            ]
        )


def main() -> None:
    args = parse_arguments()
    fractions = tuple(sorted(map(float, args.fraction.split(","))))
    args.mafft = [str(args.cores) if opt == "$CORES" else opt for opt in args.mafft]
    args.iqtree = [str(args.cores) if opt == "$CORES" else opt for opt in args.iqtree]

    logging.info("Starting process...")

    # Make output directory if not exist
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
    fract_dict = select_shared_genes(gene_dict, org_set, fractions)
    write_gene_lists(fract_dict, args.out_dir)
    align_and_trim(fract_dict, args.out_dir, args.mafft, args.trimal)

    cafiles = concat_alignments(fract_dict, args.out_dir)
    run_iqtree(cafiles, args.iqtree)

    logging.info("Job completed")


if __name__ == "__main__":
    main()
