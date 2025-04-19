#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    Busco_multigene_tree.py -i <INPUT_DIR> -f <FRACTION> [-c <CORES>]
        [-m <MAFFT_OPTION>] [-t <TRIMAL_OPTION>] [-a <AMAS_OPTION>] [-q <IQTREE_OPTION>] [-o <OUT_DIR>] [-p <PREFIX>] [-h|--help]

Purpose:
    Create multi-gene phylogenomic tree from BUSCO single-copy orthologous genes

Options:
    -h, --help                  Show this
    -i, --input_dir             Path to the input directory which contains all BUSCO outputs
    -f, --fraction              Comma-spliced fractions for creating mulit-gene phylogenetic tree
    -c, --cores                 Number of CPUs to execute [default: 8]
    -m, --mafft                 Command option for mafft [default: mafft --globalpair --maxiterate 1000 --thread $CORES]
    -t, --trimal                Command option for trimAl
    -a, --amas                  Command option for AMAS
    -q, --iqtree                Command option for IQ-TREE
    -o, --out_dir               Output directory path [default: ./]
    -p, --prefix                Prefix for output file name [default: multi-gene]

Required Packages and Softwares:
    Biopython, mafft, trimAl, AMAS, IQ-TREE

Author: Arno Hagenbeek (ArnoHagenbeek), Akito Shima (ASUQ)
Email: akito.shima@oist.jp
"""

import argparse
import logging
import subprocess
import sys
import os
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """

    # Command input using argparse
    parser = argparse.ArgumentParser(
        description="Create multi-gene phylogenomic tree from BUSCO single-copy orthologous genes.",
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
        default=[""],
        help="Command option for IQ-TREE [default: ]",
    )

    optional.add_argument(
        "-o",
        "--out_dir",
        type=Path,
        default=Path("./"),
        help="Output directory path  [default: ./]",
    )

    optional.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="multi-gene",
        help="Prefix for output file name [default: multi-gene]",
    )

    return parser.parse_args()


def validate_file(file_path: Path, valid_extensions: set) -> None:
    """
    Validates the existence and extension of a file.
    """

    if not file_path.exists():
        raise FileNotFoundError(f"File '{file_path}' does not exist.")

    if not file_path.is_file():
        raise FileNotFoundError(f"'{file_path}' is not a file.")

    if file_path.suffix.lower() not in valid_extensions:
        raise ValueError(
            f"File '{file_path}' has an invalid extension. Expected one of: {', '.join(valid_extensions)}"
        )


def run_cmd(cmd: tuple) -> int:
    """Execute a subprocess call, abort on failure

    Args:
        cmd (tuple): command

    Returns:
        int: return code
    """
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        logging.fatal(f"Error: command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(e.returncode)


def collect_gene_seqs(input_dir: Path, output_dir: Path) -> defaultdict:
    """
    Walk BUSCO outputs, aggregate single-copy faa into per-gene Fasta file.
    Returns gene_dict (gene -> [orgs]) and org_list
    """

    ## Collect sequences by gene
    gene_dict = defaultdict(set)  # gene_name: {set of orgs that have it}
    org_set = set()  # all organism names seen

    for root, _, files in os.walk(input_dir):
        for file in files:
            if not root.endswith("single_copy_busco_sequences"):
                continue

            ###
            # Need to adjust the script to collect sample ID (e.g. Genbank)
            # Maybe go upstream from the busco directory
            ###
            # Organism name is the firt subdirectory under input_dir
            rel = os.path.relpath(root, input_dir)
            org_name = rel.split(os.sep)[0]
            if org_name not in org_set:
                org_set.add(org_name)

            for fname in files:
                if not fname.endswith(".faa"):
                    continue
                gene = fname[:-4]  # strip .faa
                gene_dict[gene].add(org_name)

                in_path = os.path.join(root, fname)
                out_path = os.path.join(output_dir, f"{gene}.faa")

                with open(in_path) as inp, open(out_path, "a") as out:
                    for line in inp:
                        if line.startswith(">"):
                            out.write(f">{org_name}\n")
                        else:
                            out.write(line)

    return gene_dict, org_set


def select_shared_genes(gene_dict: dict, org_set: set, fractions: tuple) -> dict:
    """For each fraction, return list of genes present above the threshold.

    Args:
        gene_dict (dict): gene -> set(orgs)
        org_set (set): set of all organisms
        fractions (tuple): array of fractions

    Returns:
        dict: fraction -> tuple of genes present above the threshold
    """

    total = len(org_set)
    fract_dict = {}
    for frac in fractions:
        threshold = total * frac
        shared = [gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold]
        fract_dict[frac] = shared
    return fract_dict


def write_gene_lists(fract_dict: dict, output_dir: Path) -> None:
    """Write out which genes are analyzed per fraction.

    Args:
        fract_dict (dict): fraction -> tuple of genes present above the threshold
        output_dir (Path): Path of output directory
    """

    for frac, genes in fract_dict.items():
        path = os.path.join(output_dir, f"fraction{frac}_analyzed_genes.txt")
        with open(path, "w") as out:
            out.write(f"Number of genes considered: {len(genes)}\n")
            out.write("Analyzed genes:\n")
            out.write("\n".join(genes))
            out.write("\n")


def align_and_trim(
    fract_dict: dict, output_dir: Path, mafft_cmd: tuple, trimal_cmd: tuple
) -> None:
    """Align and trim all genes in the most inclusive set

    Args:
        fract_dict (dict): fraction -> tuple of genes present above the threshold
        output_dir (Path): Path of output directory
        mafft_cmd (tuple): mafft command set
        trimal_cmd (tuple): trimal command set
    """
    smallest_frac = min(fract_dict.keys())
    genes = fract_dict[smallest_frac]

    mafft_tpl = " ".join(mafft_cmd)
    trimal_tpl = " ".join(trimal_cmd)

    for gene in genes:
        infile = os.path.join(output_dir, f"{gene}.faa")
        aligned = os.path.join(output_dir, f"{gene}_aligned.faa")
        trimmed = os.path.join(output_dir, f"{gene}_trimmed.faa")
        logging.info(f"Running mafft on {infile}")
        run_cmd(mafft_tpl.split().copy() + [infile, ">", aligned])
        logging.info(f"Running trimAl on {infile}")
        run_cmd(trimal_tpl.split().copy() + ["-in", aligned, "-out", trimmed])


def concat_and_build(fract_dict: dict, output_dir: Path, iqtree_cmd: tuple) -> None:
    """For each fraction, concatenate trimmed gene alignments and run IQ-TREE

    Args:
        fract_dict (dict): fraction -> tuple of genes present above the threshold
        output_dir (Path): Path of output directory
        iqtree_cmd (tuple): IQ-TREE command set
    """
    iq_tpl = " ".join(iqtree_cmd)

    for frac, genes in fract_dict.items():
        concat = defaultdict(dict)
        for gene in genes:
            path = os.path.join(output_dir, f"{gene}_trimmed.faa")
            org = None
            seq = []
            with open(path) as inp:
                for line in inp:
                    if line.startswith(">"):
                        org = line[1:].strip()
                    else:
                        seq.append(line.strip())
            concat[org][gene] = "".join(seq)

        # pad missing genes
        for org in concat:
            for gene in genes:
                if gene not in concat[org]:
                    length = len(next(iter(concat.values()))[gene])
                    concat[org][gene] = "-" * length

        # Write concatenated file
        cafile = os.path.join(output_dir, f"fraction{frac}_concat.faa")
        with open(cafile, "w") as out:
            for org, gdict in concat.items():
                out.write(f">{org}\n")
                out.write("".join(gdict[g] for g in genes) + "\n")

        # Lastly, run iqtree on the concatenated file (per given fraction)
        logging.info(f"Running IQ-TREE on {cafile}")
        cmd = iq_tpl.split().copy() + [
            "-s",
            cafile,
            "-pre",
            f"fraction{frac}_concat.faa",
        ]
        run_cmd(cmd)


def main() -> None:
    ## Parse arguments
    args = parse_arguments()
    fractions = tuple(sorted(list(map(float, args.fraction.split(",")))))

    logging.info("Starting process...")

    # Make output directory if not exist
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
    fract_dict = select_shared_genes(gene_dict, org_set, fractions)
    write_gene_lists(fract_dict, args.out_dir)
    align_and_trim(fract_dict, args.out_dir, args.mafft, args.trimal)
    concat_and_build(fract_dict, args.out_dir, args.iqtree)


if __name__ == "__main__":
    main()
