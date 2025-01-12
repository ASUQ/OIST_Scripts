#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    hyp_protein_extractor -i <INPUT> [-o <OUT>] [-p <PREFIX>] [-h|--help]

Purpose:
	Extract hypothetical proteins from faa file (Prokka annotated)

Options:
    -h, --help                  Show this
    -i, --input                 faa file of annotated protein - generated from prokka
    -o, --out <OUT>             Output directory path (string) [default: ./]
    -p, --prefix <PREFIX>       Prefix for output file. [default: hyp_protein]

Required Packages:
    Biopython

Author: Akito Shima (ASUQ)
Email: akito-shima@oist.jp
"""

def parse_arguments():
    """
    Parse command-line arguments.
    """

    import argparse

    # Command input using argparse
    parser = argparse.ArgumentParser(
        description= 'Extract hypothetical proteins from faa file (Prokka annotated)',
        epilog='Required package: Biopython',
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-i", "--input",
		type=str,
		required=True,
		help="faa file of annotated protein - generated from prokka"
	)

    optional.add_argument(
        "-o", "--out",
        type=str,
        default='./',
        help="Output directory path [default: ./]"
    )

    optional.add_argument(
        "-p", "--prefix",
        type=str,
        default='hyp_protein',
        help="Prefix for output file [default: hyp_protein]"
    )

    return parser.parse_args()


def check_required_packages(packages):
    """
    Check if the required packages are installed.
    """
    import importlib.util
    import sys

    missing_packages = []
    for package in packages:
        if importlib.util.find_spec(package) is None:
            missing_packages.append(package)

    if missing_packages:
        print(f"Error: The following required packages are missing: {', '.join(missing_packages)}")
        sys.exit(1)


def validate_file(file_path, valid_extensions):
    """
    Validates the existence and extension of a file.
    """
    import os

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' does not exist.")

    root, ext = os.path.splitext(file_path)
    if ext.lower() not in valid_extensions:
        raise ValueError(f"File '{file_path}' has an invalid extension. Expected one of: {', '.join(valid_extensions)}")


def hyp_extractor(gene, output_path):
    from Bio import SeqIO

    # Extract sequences of the targets and save as fasta file
    sequence_records = list()
    for record in SeqIO.parse(gene, 'fasta'):
        if 'hypothetical protein' in record.description:
            sequence_records.append(record)

    if len(sequence_records):
        SeqIO.write(sequence_records, output_path, 'fasta')
    else:
        print('No hypothetical protein sequence was found.')


def main():
    # Parse arguments
    args = parse_arguments()

	# Check the installed packages
    required_packages = ["Bio"]
    check_required_packages(required_packages)

	# Package import
    from Bio import SeqIO
    import os.path
    import sys

	# Extract argument values
    gene = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix

	# Check file existence and extension
    validate_file(gene, ['.fasta', '.fas', '.fa', '.faa', '.mpfa'])

    # Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

	# Set the path for the output file
    output_file_path = os.path.join(output_dir, f"{prefix}.faa")
    if os.path.exists(output_file_path):
        response = input(f"Output file '{output_file_path}' already exists. Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("Aborting to avoid overwriting the existing file.")
            sys.exit(1)

    hyp_extractor(gene, output_file_path)


if __name__ == '__main__':
    main()
