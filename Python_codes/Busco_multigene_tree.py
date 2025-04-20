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
    -o, --out_dir               Output directory path [default: ./output]

Required Packages and Softwares:
    Biopython, mafft, trimAl, AMAS, IQ-TREE

Author: Arno Hagenbeek (ArnoHagenbeek), Akito Shima (ASUQ)
Email: akito.shima@oist.jp
"""

import argparse
import logging
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

SOFTWARES = ["mafft", "trimal", "AMAS.py", "iqtree"]


def check_software(tool: str) -> None:
    """Verify that a required software is available in PATH."""
    if shutil.which(tool) is None:
        logging.fatal(f"Required software '{tool}' not found in PATH.")
        sys.exit(1)


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
        help="Command option for mafft [default: mafft --globalpair --maxiterate 1000 --thread $CORES $INPUT > $OUTPUT]",
    )

    optional.add_argument(
        "-t",
        "--trimal",
        nargs="+",
        metavar="TRIMAL_OPTION",
        default=["--automated1"],
        help="Command option for trimal [default: trimal --automated1 -in $INPUT -out $OUTPUT]",
    )

    optional.add_argument(
        "-a",
        "--amas",
        nargs="+",
        metavar="AMAS_OPTION",
        default=[
            "concat",
            "--in-format",
            "fasta",
            "--cores",
            "$CORES",
            "--data-type",
            "aa",
            "--part-format",
            "nexus",
        ],
        help="Command option for AMAS [default: AMAS.py concat --in-files $INPUT --in-format fasta --data-type aa --concat-out $OUT_concat.faa --concat-part $OUT_partitions.txt --part-format nexus --cores $Cores]",
    )

    optional.add_argument(
        "-q",
        "--iqtree",
        nargs="+",
        metavar="IQTREE_OPTION",
        default=["-B", "1000", "-alrt", "1000", "-T", "$CORES"],
        help="Command option for IQ-TREE [default: iqtree -B 1000 -alrt 1000 -T $CORES -s $INPUT -p $PARTITION -pre $PREFIX]",
    )

    optional.add_argument(
        "-o",
        "--out_dir",
        type=Path,
        default=Path("./output"),
        help="Output directory path  [default: ./output]",
    )

    return parser.parse_args()


def run_cmd(cmd: list[str], stdout=None) -> None:
    """Execute a subprocess call, abort on failure"""
    try:
        subprocess.check_call(cmd, stdout=stdout)
    except subprocess.CalledProcessError as e:
        logging.fatal(f"Command failed {e.returncode}: {' '.join(cmd)}")
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
        # assuming structure: /input_dir/subdirectories/{sample_name}/busco_output/run_lineage/single_copy_busco_sequences/*.faa
        ###
        try:
            org_name = seq_dir.relative_to(input_dir).parts[-4]
        except IndexError:
            raise ValueError(
                f"Directory {seq_dir} too shallow for extracting organism name."
            )

        if org_name in org_set:
            raise ValueError(
                f"Organism Name is duplicated. Please check the directory structure."
            )
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
    """For each fraction, return list of genes present above the threshold."""

    total = len(org_set)
    frac_dict: dict[float, list[str]] = {}
    for frac in fractions:
        threshold = math.ceil(total * frac)
        frac_dict[frac] = [
            gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold
        ]
    return frac_dict


def write_gene_lists(frac_dict: dict[float, list[str]], output_dir: Path) -> None:
    """Dump which genes were used at each completeness threshold."""

    for frac, genes in frac_dict.items():
        file_path = output_dir / f"frac{int(frac * 100)}pct_genes.txt"
        with file_path.open("w") as out:
            out.write(f"Number of genes considered: {len(genes)}\n")
            out.write("Analyzed genes:\n")
            out.write("\n".join(genes) + "\n")


def align_and_trim(
    frac_dict: dict[float, list[str]],
    output_dir: Path,
    mafft_opts: list[str],
    trimal_opts: list[str],
) -> None:
    """Align and trim all genes in the most inclusive set"""

    smallest_frac = min(frac_dict.keys())
    genes = frac_dict[smallest_frac]

    for gene in genes:
        infile = output_dir / f"{gene}.faa"
        aligned = output_dir / f"{gene}_aligned.faa"
        trimmed = output_dir / f"{gene}_trimmed.faa"

        logging.info(f"Running mafft {str(infile)} -> {str(aligned)}")
        with aligned.open("w") as out_f:
            run_cmd(["mafft"] + mafft_opts + [str(infile)], stdout=out_f)

        logging.info(f"Running trimAl {str(aligned)} -> {str(trimmed)}")
        run_cmd(["trimal"] + trimal_opts + ["-in", str(aligned), "-out", str(trimmed)])


def concat_alignments(
    frac_dict: dict[float, list[str]], output_dir: Path, amas_opts: list[str]
) -> dict[float, tuple[Path, Path]]:
    """Run AMAS to concatenate trimmed gene alignments using AMAS"""

    cafiles: dict[float, tuple[Path, Path]] = {}
    for frac, genes in frac_dict.items():
        if not genes:
            logging.warning(f"No genes for fraction {frac}, skipping...")
            continue

        trimmed_files = [str(output_dir / f"{gene}_trimmed.faa") for gene in genes]
        concat_faa = output_dir / f"frac{int(frac * 100)}pct_concat.faa"
        partition_file = output_dir / f"frac{int(frac * 100)}pct_partitions.nex"

        logging.info(f"Running AMAS concat: fraction {frac}")

        cmd = (
            ["AMAS.py"]
            + amas_opts
            + [
                "--concat-out",
                str(concat_faa),
                "--concat-part",
                str(partition_file),
                "--in-files",
            ]
            + trimmed_files
        )

        run_cmd(cmd)
        cafiles[frac] = (concat_faa, partition_file)

    return cafiles


def run_iqtree(
    cafiles: dict[float, tuple[Path, Path]], output_dir: Path, iqtree_opts: list[str]
) -> None:
    """Run IQ-TREE"""

    for frac, cafile in cafiles.items():
        concat_faa, partition_file = cafile

        results_dir = output_dir / f"frac{int(frac * 100)}pct_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        tree_prefix = results_dir / f"frac{int(frac * 100)}pct"

        logging.info(f"Running IQ-TREE on {str(concat_faa)} and {str(partition_file)}")
        run_cmd(
            ["iqtree"]
            + iqtree_opts
            + [
                "-s",
                str(concat_faa),
                "-p",
                str(partition_file),
                "-pre",
                str(tree_prefix),
            ]
        )


def main() -> None:
    # Verify software in the PATH
    for tool in SOFTWARES:
        check_software(tool)

    args = parse_arguments()

    try:
        args.out_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        logging.error(f"Output directory exists")
        sys.exit(1)

    # Parse fractions, ensure 0<frac<=1
    fractions = set()
    for frac in args.fraction.split(","):
        try:
            frac = float(frac)
        except ValueError:
            logging.error("fractions must be comma-spliced numbers")
            sys.exit(1)
        if not 0 < frac <= 1:
            logging.error("fractions must be numbers between 0 and 1")
            sys.exit(1)

        fractions.add(round(frac, 2))
    fractions = tuple(sorted(fractions))

    # Replacing $CORES with exact CPU number selected by the user
    args.mafft = [opt.replace("$CORES", str(args.cores)) for opt in args.mafft]
    args.amas = [opt.replace("$CORES", str(args.cores)) for opt in args.amas]
    args.iqtree = [opt.replace("$CORES", str(args.cores)) for opt in args.iqtree]

    logging.info("Starting process...")

    gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
    frac_dict = select_shared_genes(gene_dict, org_set, fractions)
    write_gene_lists(frac_dict, args.out_dir)
    align_and_trim(frac_dict, args.out_dir, args.mafft, args.trimal)

    cafiles = concat_alignments(frac_dict, args.out_dir, args.amas)
    run_iqtree(cafiles, args.out_dir, args.iqtree)

    logging.info("Job completed")


if __name__ == "__main__":
    main()
