#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    fasta_to_tsv -i <INPUT> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Purpose:
	Convert FASTA file records to a TSV file (name, length, coverage)
	This code adds dummy coverage (1) for the use in blobtools

Options:
    -h, --help                  Show this
    -i, --input <INPUT>         Fasta file
    -p, --prefix <PREFIX>       Prefix for output file.
    -o, --out <OUT>             Output directory path [default: ./]

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
        description='Convert FASTA file records to a TSV file (name, length, coverage).',
        epilog='required package: Biopython',
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-i", "--input",
		type=str,
		required=True,
		help="Fasta file"
    )

    optional.add_argument(
        "-o", "--out",
		type=str,
		default='./',
		help='Output directory path  [default: ./]'
    )

    optional.add_argument(
        "-p", "--prefix",
		type=str,
		help='Prefix for output file'
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


def fasta_to_tsv(input_path, output_path):
    from Bio import SeqIO
    import os

    records = SeqIO.parse(input_path, 'fasta')

    record_list = []

    for record in records:
        name = record.id
        length = len(record.seq)
        cov = 1
        record_list.append((name, length, cov))

    with open(output_path, 'w') as output_handle:
        output_handle.write(f'#name\tlength\tcoverage\n')
        for record in record_list:
            output_handle.write(f"{record.id}\t{len(record.seq)}\t1\n")


def main():
    # Parse arguments
    args = parse_arguments()

	# Check the installed packages
    required_packages = ['Bio']
    check_required_packages(required_packages)

	# Package import
    from Bio import SeqIO
    import os
    import sys

	# Extract argument values
    input_path = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix

	# Check file existence and extension
    validate_file(input_path, ['.fasta', '.fas', '.fa', '.fna', '.mpfa'])

    # Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Check prefix name
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(input_path))[0]

	# Set the path for the output file
    output_file_path = os.path.join(output_dir, f"{prefix}.tsv")
    if os.path.exists(output_file_path):
        response = input(f"Output file '{output_file_path}' already exists. Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print('Aborting to avoid overwriting the existing file.')
            sys.exit(1)

	#==========================================================#

	# Convert fasta records to TSV
    fasta_to_tsv(input_path, output_file_path)

if __name__ == '__main__':
	main()
