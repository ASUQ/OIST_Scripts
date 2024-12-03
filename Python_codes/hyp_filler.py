#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    hyp_filler -a <ANNOTATION> -l <HYP> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Purpose:
	??? (Forgot why I wrote this...)

Options:
    -h, --help                  Show this
    -a, --annotation            faa file of annotated protein - generated from prokka
    -l, --hyp                   List of Loci which are hyp
    -p, --prefix <PREFIX>       Prefix for output file name [default: hyp_filler]
    -o, --out <OUT>             Output directory path (string) [default: ./]

Required Packages:
    Biopython, pandas

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
        epilog="Required package: Biopython, pandas",
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-a", "--annotation",
		type=str,
		required=True,
		help="faa file of annotated protein - generated from prokka"
	)

    required.add_argument(
        "-l", "--hyp",
		type=str,
		required=True,
		help="List of Loci which are hyp"
	)

    optional.add_argument(
        "-o", "--out",
		type=str,
		default='./',
		help="Output directory path  [default: ./]"
	)

    optional.add_argument(
        "-p", "--prefix",
		type=str,
		default='hyp_filled',
		help="Prefix for output file name [default: hyp_filled]"
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


def hyp_filler(annotation, hyp, output_path):
    from Bio import SeqIO
    import pandas as pd

    # Rename the hypothetical proteins
    hyp_df = pd.read_csv(hyp)
    hyp_id = list(hyp_df['LOCUS'])

    sequence_records = list()
    for record in SeqIO.parse(annotation, 'fasta'):
        if record.id in hyp_id:
            new_record = record
            new_record.description = record.id + '_' + str(list(hyp_df[hyp_df['LOCUS'] == record.id]['protein'])[0])
            record = new_record
        sequence_records.append(record)

    SeqIO.write(sequence_records, output_path, 'fasta')


def main():
    # Parse arguments
    args = parse_arguments()

	# Check the installed packages
    required_packages = ['Bio', 'pandas']
    check_required_packages(required_packages)

	# Package import
    from Bio import SeqIO
    import os.path
    import pandas as pd
    import sys

    # Parse the inputs from command-line
    annotation = os.path.abspath(args.annotation)
    hyp = os.path.abspath(args.hyp)
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix

	# Check file existence and extension
    validate_file(annotation, ['.fasta', '.fas', '.fa', '.faa', '.mpfa'])
    validate_file(hyp, ['.csv'])

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


    hyp_filler(annotation, hyp, output_file_path)

if __name__ == '__main__':
    main()
