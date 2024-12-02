#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    sequence_extractor -i <INPUT> -f <FASTA> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Purpose:
	Extract records from a FASTA file based on a list of record_id

Options:
    -h, --help                  Show this
    -i, --input                 List of seqeunces you want (string) - list by "\n"
    -f, --fasta <FASTA>         Fasta file of sequences (string)
    -p, --prefix <PREFIX>       Prefix for output file.
    -o, --out <OUT>             Output directory path (string) [default: ./]

Required Packages:
    Biopython

Author: Akito Shima (ASUQ)
Email: akito-shima@oist.jp
"""

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
        print("Please install them using the following command:")
        print(f"    conda install -c conda-forge {' '.join(missing_packages)}")
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


def find_duplicates(seq_list):
    """
    Finds duplicate items in a list.

    Parameters:
        seq_list (list): List of sequence IDs.

    Returns:
        duplicates (set): Set of duplicate sequence IDs.
    """
    seen = set()
    duplicates = set()
    for seq_id in seq_list:
        if seq_id in seen:
            duplicates.add(seq_id)
        else:
            seen.add(seq_id)
    return duplicates


def main():
	# Check the installed packages
    required_packages = ["Bio"]
    check_required_packages(required_packages)

	# Package import
    import argparse
    from Bio import SeqIO
    import os.path
    import sys

    # Command input using argparse
    parser = argparse.ArgumentParser(
        description="Extract records from a FASTA file based on a list of record_id.",
        epilog="required package: Biopython"
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="List of sequences you want (string)"
    )

    required.add_argument(
        "-f", "--fasta",
        type=str,
        required=True,
        help="Fasta file of sequences (string)"
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
        help="Prefix for output file  [default: input file name]"
    )

    # Parse the inputs from command-line
    args = parser.parse_args()
    nodes = args.input
    sequences = args.fasta
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix


	# Check file existence and extension
    validate_file(nodes, ['.txt'])
    validate_file(sequences, ['.fasta', '.fas', '.fa', '.fna'])

    # Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Check prefix name
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(nodes))[0]

	# Set the path for the output file
    output_file_path = os.path.join(output_dir, f"{prefix}.fasta")
    if os.path.exists(output_file_path):
        response = input(f"Output file '{output_file_path}' already exists. Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("Aborting to avoid overwriting the existing file.")
            sys.exit(1)

	#==========================================================#

	# Save the sequence ids in names list
    with open(nodes, 'r') as f:
        names = []
        for line in f:
            seq_id = line.strip()
            if seq_id:
                names.append(seq_id)

        duplicates = find_duplicates(names)
        if duplicates:
            print(f"Warning: {len(duplicates)} duplicate sequence IDs found in the input file: {', '.join(duplicates)}")

        if not names:
            print("Input file is empty. Please provide a file with sequence IDs.")
            sys.exit(1)


    # Extract sequences of the targets and save as fasta file
    output_handle = None
    found_ids = set()
    with open(sequences, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in names:
                if not output_handle:
                    output_handle = open(output_file_path, 'w')

                SeqIO.write(record, output_handle, 'fasta')
                found_ids.add(record.id)

        if output_handle:
            output_handle.close()


    unmatched_ids = [name for name in names if name not in found_ids]
    if unmatched_ids:
        print(f"Warning: The following sequence IDs were not found in the FASTA file: {', '.join(unmatched_ids)}")

    if output_handle:
        print(f"Sequences have been written to {output_file_path}")

    else:
        print('No sequence were found.')
        sys.exit(1)

if __name__ == '__main__':
	main()