#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    Busco_multigene_tree.py <subcommand> [options]

Purpose:
    Create multi-gene phylogenomic tree from BUSCO single-copy orthologs.

Subcommands:
    all                         Run all steps: collect, select, align, infer
    collect                     Collect per-gene FASTA files from BUSCO outputs
    select                      Select shared genes & write gene lists
    align                       Align & trim gene alignments (mafft & trimAl)
    infer                       Concatenate & build tree (AMAS & IQ-TREE)

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
    Biopython, tqdm, mafft, trimAl, AMAS, IQ-TREE

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
from tqdm import tqdm


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

SOFTWARES = ["mafft", "trimal", "AMAS.py", "iqtree"]

OPTION_DEFS = {
    "input_dir": dict(
        flags=("-i", "--input_dir"),
        kwargs=dict(type=Path, required=True, help="Path to BUSCO outputs"),
    ),
    "out_dir": dict(
        flags=("-o", "--out_dir"),
        kwargs=dict(
            type=Path,
            default=Path("./output"),
            help="Output directory [default: ./output]",
        ),
    ),
    "cores": dict(
        flags=("-c", "--cores"),
        kwargs=dict(type=int, default=8, help="CPUs to use [default: 8]"),
    ),
    "verbose": dict(
        flags=("-v", "--verbose"),
        kwargs=dict(action="store_true", help=argparse.SUPPRESS),
    ),
    "fraction": dict(
        flags=("-f", "--fraction"),
        kwargs=dict(
            type=str,
            default="0.9",
            help="Comma-spliced completeness fraction(s) (e.g. '0.8,0.9') [default: 0.9]",
        ),
    ),
    "mafft": dict(
        flags=("-m", "--mafft"),
        kwargs=dict(
            nargs="+",
            metavar="MAFFT_OPTION",
            default=["--globalpair", "--maxiterate", "1000", "--thread", "$CORES"],
            help="MAFFT options [default: mafft --globalpair --maxiterate 1000 --thread $CORES $INPUT > $OUTPUT]",
        ),
    ),
    "trimal": dict(
        flags=("-t", "--trimal"),
        kwargs=dict(
            nargs="+",
            metavar="TRIMAL_OPTION",
            default=["--automated1"],
            help="trimAl options [default: trimal --automated1 -in $INPUT -out $OUTPUT]",
        ),
    ),
    "amas": dict(
        flags=("-a", "--amas"),
        kwargs=dict(
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
            help="AMAS options [default: AMAS.py concat --in-files $INPUT --in-format fasta --data-type aa --concat-out $OUT_concat.faa --concat-part $OUT_partitions.txt --part-format nexus --cores $Cores]",
        ),
    ),
    "iqtree": dict(
        flags=("-q", "--iqtree"),
        kwargs=dict(
            nargs="+",
            metavar="IQTREE_OPTION",
            default=["-B", "1000", "-alrt", "1000", "-T", "$CORES"],
            help="IQ-TREE options [default: iqtree -B 1000 -alrt 1000 -T $CORES -s $INPUT -p $PARTITION -pre $PREFIX]",
        ),
    ),
}


SUBCMD_OPTS = {
    "collect": (
        ("input_dir", "out_dir", "verbose"),
        "Collect per-gene FASTA files from BUSCO outputs",
    ),
    "select": (
        ("input_dir", "out_dir", "fraction", "verbose"),
        "Select shared genes & write gene lists",
    ),
    "align": (
        ("input_dir", "out_dir", "fraction", "cores", "mafft", "trimal", "verbose"),
        "Align & trim gene alignments (mafft & trimAl)",
    ),
    "infer": (
        ("input_dir", "out_dir", "fraction", "cores", "amas", "iqtree", "verbose"),
        "Concatenate & build tree (AMAS & IQ-TREE)",
    ),
    "all": (
        (
            "input_dir",
            "out_dir",
            "fraction",
            "cores",
            "mafft",
            "trimal",
            "amas",
            "iqtree",
            "verbose",
        ),
        "Run all steps: collect, select, align, infer",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build multi-gene phylogenomic tree from BUSCO single-copy orthologous gene",
        epilog="Required package and Softwares: Biopython, mafft, trimAl, AMAS, IQ-TREE",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subs = parser.add_subparsers(dest="command", required=True)

    for cmd, (opt_keys, help) in SUBCMD_OPTS.items():
        p = subs.add_parser(cmd, help=f"{help}")
        for key in opt_keys:
            flags = OPTION_DEFS[key]["flags"]
            kwargs = OPTION_DEFS[key]["kwargs"]
            p.add_argument(*flags, **kwargs)

    return parser.parse_args()


def check_software(tool: str) -> None:
    """Verify that a required software is available in PATH."""
    if shutil.which(tool) is None:
        logging.fatal(f"Required software '{tool}' not found in PATH.")
        sys.exit(1)


def run_cmd(cmd: list[str], stdout=None) -> None:
    """Execute a subprocess call, abort on failure"""
    try:
        logging.debug(f"Executing: {cmd}")
        subprocess.check_call(cmd, stdout=stdout)
    except subprocess.CalledProcessError as e:
        logging.fatal(f"Command failed {e.returncode}: {' '.join(cmd)}")
        sys.exit(e.returncode)


## Functions for parsing arguments
def parse_fractions(frac_str: str) -> list[float]:
    """Parse comma-separated fractions into a sorted list of floats."""
    fracs = set()
    for frac in frac_str.split(","):
        try:
            frac = float(frac)
        except ValueError:
            logging.error("fractions must be comma-spliced numbers")
            sys.exit(1)
        if not 0 < frac <= 1:
            logging.error("fractions must be numbers between 0 and 1")
            sys.exit(1)

        fracs.add(round(frac, 2))

    return sorted(fracs)


def load_genes_for_fraction(frac: float, output_dir: Path) -> list[str]:
    """Load the list of selected genes for a given fraction from disk."""
    pct = int(frac * 100)
    path = output_dir / f"frac{pct}pct_results/frac{pct}pct_genes.txt"
    with path.open() as f:
        lines = f.read().splitlines()
    return lines[2:]  # skip header lines


def collect_gene_seqs(
    input_dir: Path, output_dir: Path
) -> tuple[dict[str, set[str]], set[str]]:
    """
    Parse single-copy BUSCO FASTAs with SeqIO, prefix record.id with organism name,
    and append into per-gene files.
    """

    logging.info("Collecting genes from BUSCO outputs")

    raw_dir = output_dir / "seqs" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    # gene_name: {set of orgs that have it}
    gene_dict: dict[str, set[str]] = defaultdict(set)
    org_set: set[str] = set()  # all organism names seen

    for seq_dir in input_dir.rglob(
        "single_copy_busco_sequences", recurse_symlinks=True
    ):
        logging.debug(f"Collecting genes from {seq_dir}")

        if not seq_dir.is_dir():
            logging.debug(f"{seq_dir} is not directory")
            continue

        ###
        # Extract sample name
        # assuming structure: /input_dir/subdirectories/{sample_name}/busco_output/run_lineage/busco_sequences/single_copy_busco_sequences/*.faa
        ###
        try:
            org_name = seq_dir.relative_to(input_dir).parts[-5]
        except IndexError:
            raise ValueError(
                f"Directory {seq_dir} too shallow for extracting organism name."
            )

        if org_name in org_set:
            raise ValueError(
                f"Organism Name is duplicated. Please check the directory structure."
            )

        logging.debug(f"org_name: {org_name} found")
        org_set.add(org_name)

        logging.debug(f"Extracting genes from faa files")
        for faa_file in seq_dir.glob("*.faa"):
            logging.debug(f"Extracting gene seq from {str(faa_file)}")

            gene = faa_file.stem
            gene_dict[gene].add(org_name)

            out_file = raw_dir / f"{gene}.faa"
            with faa_file.open("r") as inp, out_file.open("a") as out:
                for rec in SeqIO.parse(inp, "fasta"):
                    rec.id = org_name
                    rec.description = ""
                    SeqIO.write(rec, out_file, "fasta")

    logging.debug(f"org_set: {org_set}")
    logging.debug(f"gene_dict: {gene_dict}")

    return gene_dict, org_set


def select_shared_genes(
    gene_dict: dict[str, set[str]], org_set: set[str], fractions: list[float]
) -> dict[float, list[str]]:
    """For each fraction, return list of genes present above the threshold."""

    logging.info("Selecting shared genes")

    total = len(org_set)
    frac_dict: dict[float, list[str]] = {}
    for frac in fractions:
        threshold = math.ceil(total * frac)
        frac_dict[frac] = [
            gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold
        ]

    logging.debug(f"frac_dict: {frac_dict}")
    return frac_dict


def write_gene_lists(frac_dict: dict[float, list[str]], output_dir: Path) -> None:
    """Write out which genes pass completeness threshold."""

    logging.info("Writing out gene lists")

    for frac, genes in frac_dict.items():
        pct = int(frac * 100)
        results_dir = output_dir / f"frac{pct}pct_results"
        results_dir.mkdir(parents=True, exist_ok=True)

        file_path = results_dir / f"frac{pct}pct_genes.txt"
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

    logging.info("Aligning the genes")

    seq_dir = output_dir / "seqs"
    raw_dir = seq_dir / "raw"
    aligned_dir = seq_dir / "aligned"
    trimmed_dir = seq_dir / "trimmed"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    smallest_frac = min(frac_dict.keys())
    genes = frac_dict[smallest_frac]

    logging.debug(
        f"Aligning genes in fraction: {smallest_frac} which has {len(genes)} genes"
    )

    for gene in tqdm(genes, desc="Align & trim"):
        infile = raw_dir / f"{gene}.faa"
        aligned = aligned_dir / f"{gene}_aligned.faa"
        trimmed = trimmed_dir / f"{gene}_trimmed.faa"

        logging.info(f"Running mafft {str(infile)} -> {str(aligned)}")
        with aligned.open("w") as out_f:
            run_cmd(["mafft"] + mafft_opts + [str(infile)], stdout=out_f)

        logging.info(f"Running trimAl {str(aligned)} -> {str(trimmed)}")
        run_cmd(["trimal"] + trimal_opts + ["-in", str(aligned), "-out", str(trimmed)])


def concat_alignments(
    frac_dict: dict[float, list[str]], output_dir: Path, amas_opts: list[str]
) -> dict[float, tuple[Path, Path]]:
    """Run AMAS to concatenate trimmed gene alignments using AMAS"""

    logging.info(f"Concatenating alignments")

    cafiles: dict[float, tuple[Path, Path]] = {}
    for frac, genes in frac_dict.items():
        if not genes:
            logging.warning(f"No genes for fraction {frac}, skipping...")
            continue

        pct = int(frac * 100)
        results_dir = output_dir / f"frac{pct}pct_results"
        concat_faa = results_dir / "concat.faa"
        partition_file = results_dir / "partitions.nex"

        trimmed_files = [
            str(output_dir / f"seqs/trimmed/{g}_trimmed.faa") for g in genes
        ]

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

    logging.info("Building trees")

    for frac, (concat_faa, partition_file) in cafiles.items():
        pct = int(frac * 100)

        prefix = output_dir / f"frac{pct}pct_results/frac{pct}pct"

        logging.info(f"Running IQ-TREE on {str(concat_faa)} and {str(partition_file)}")
        run_cmd(
            ["iqtree"]
            + iqtree_opts
            + ["-s", str(concat_faa), "-p", str(partition_file), "-pre", str(prefix)]
        )


def main() -> None:
    # Verify software in the PATH
    for tool in SOFTWARES:
        check_software(tool)

    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logging.debug(f"parsed arguments: {args!r}")

    if not args.input_dir.is_dir():
        logging.error(f"Input {str(args.input_dir)} is not directory")
        sys.exit(1)

    logging.info("Starting process...")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.command == "collect":
        logging.info("Running subcomannd: collect")
        collect_gene_seqs(args.input_dir, args.out_dir)

    elif args.command == "select":
        logging.info("Running subcomannd: select")
        fracs = parse_fractions(args.fraction)
        gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
        frac_dict = select_shared_genes(gene_dict, org_set, fracs)
        write_gene_lists(frac_dict, args.out_dir)

    elif args.command == "align":
        logging.info("Running subcomannd: align")
        fracs = parse_fractions(args.fraction)
        args.mafft = [opt.replace("$CORES", str(args.cores)) for opt in args.mafft]

        gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
        frac_dict = select_shared_genes(gene_dict, org_set, fracs)
        align_and_trim(frac_dict, args.out_dir, args.mafft, args.trimal)

    elif args.command == "infer":
        logging.info("Running subcomannd: infer")
        fracs = parse_fractions(args.fraction)
        args.amas = [opt.replace("$CORES", str(args.cores)) for opt in args.amas]
        args.iqtree = [opt.replace("$CORES", str(args.cores)) for opt in args.iqtree]

        frac_dict = {}
        for frac in fracs:
            genes = load_genes_for_fraction(frac, args.out_dir)
            frac_dict[frac] = genes
        cafiles = concat_alignments(frac_dict, args.out_dir, args.amas)
        run_iqtree(cafiles, args.out_dir, args.iqtree)

    elif args.command == "all":
        logging.info("Running subcomannd: all")
        fracs = parse_fractions(args.fraction)
        args.mafft = [opt.replace("$CORES", str(args.cores)) for opt in args.mafft]
        args.amas = [opt.replace("$CORES", str(args.cores)) for opt in args.amas]
        args.iqtree = [opt.replace("$CORES", str(args.cores)) for opt in args.iqtree]

        gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir)
        frac_dict = select_shared_genes(gene_dict, org_set, fracs)
        write_gene_lists(frac_dict, args.out_dir)
        align_and_trim(frac_dict, args.out_dir, args.mafft, args.trimal)
        cafiles = concat_alignments(frac_dict, args.out_dir, args.amas)
        run_iqtree(cafiles, args.out_dir, args.iqtree)

    else:
        logging.error(f"Unknown command {args.command}")
        sys.exit(1)

    logging.info("Job completed")


if __name__ == "__main__":
    main()
