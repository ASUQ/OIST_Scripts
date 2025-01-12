#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    hyp_protein_table -i1 <BLASTP> -i2 <PHMMER> [-o <OUT>] [-p <PREFIX>] [-h|--help]

Purpose:
	Create the table which list the search result from blastp and phmmer

Options:
    -h, --help                  Show this
    -i1, --input1 <BLASTP>      tsv file from blastp
    -i2, --input2 <PHMMER>      tsv file from phmmer
    -o, --out <OUT>             Output directory path (string) [default: ./]
    -p, --prefix <PREFIX>       Prefix for output file. [default: hyp_protein_table]

Required Packages:
    pandas

Author: Akito Shima (ASUQ)
Email: akito-shima@oist.jp
"""

MAX_PROTEINS = 5

def parse_arguments():
    """
    Parse command-line arguments.
    """

    import argparse

    # Command input using argparse
    parser = argparse.ArgumentParser(
        description= 'Create the table which list the search result from blastp and phmmer',
        epilog='Required package: Biopython, pandas',
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-i1", "--input1",
		type=str,
		required=True,
		help='tsv file from blastp'
	)

    required.add_argument(
        "-i2", "--input2",
		type=str,
		required=True,
		help='tsv file from phmmer'
    )

    optional.add_argument(
        "-o", "--out",
		type=str,
		default='./',
		help='Output directory path [default: ./]'
	)

    optional.add_argument(
        "-p", "--prefix",
		type=str,
		default='hyp_protein_table',
		help='Output directory path [default: hyp_protein_table]'
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


def load_and_extract_proteins(file_path, usecols, column_name, regex_pattern):
	"""
    Load a tsv file and extract protein names using a regex pattern
	"""
	df = pd.read_csv(file_path, usecols=usecols)
	df['protein_name'] = df[column_name].str.extract(regex_pattern)
	return df


def collect_protein_results(df, locus_col, protein_col):
    """
    Collect protein names per locus without duplication and limit to MAX_PROTEINS.
    """
    protein_dict = defaultdict(list)
    for locus, protein in zip(df[locus_col], df[protein_col]):
        if pd.notna(protein) and protein not in protein_dict[locus]:
            protein_dict[locus].append(protein)

    # Ensure each list has exactly MAX_PROTEINS items
    for locus in protein_dict:
        proteins = protein_dict[locus]
        proteins = (proteins + ['nan'] * MAX_PROTEINS)[:MAX_PROTEINS]
        protein_dict[locus] = ';'.join(proteins)

    # Convert to DataFrame
    result_df = pd.DataFrame.from_dict(protein_dict, orient='index', columns=[protein_col])
    result_df.index.name = 'LOCUS'
    return result_df


def write_results(result_df, result_path, final_path):
    """
    Write the results to CSV files as per the original script's output format.
    """
    # Write the first result file
    with open(result_path, mode='w', newline='\n') as f:
        f.write('LOCUS,blastp,phmmer\n')
        for locus, row in result_df.iterrows():
            blastp_proteins = row['blastp_protein'].split(';')
            phmmer_proteins = row['phmmer_protein'].split(';')
            for bp, ph in zip(blastp_proteins, phmmer_proteins):
                f.write(f'{locus},{bp},{ph}\n')
            f.write('\n')

    # Write the final result file
    with open(final_path, mode='w', newline='\n') as f:
        f.write('LOCUS,protein,gene\n')
        for locus in result_df.index:
            f.write(f'{locus},,\n')


def main():
    # Parse arguments
    args = parse_arguments()

	# Check the installed packages
    required_packages = ['pandas']
    check_required_packages(required_packages)

	# Package import
    from Bio import SeqIO
    from collections import defaultdict
    import os.path
    import pandas as pd
    import re
    import sys

	# Extract argument values
    blastp = os.path.abspath(args.input1)
    phmmer = os.path.abspath(args.input2)
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix

	# Check file existence and extension
    validate_file(blastp, ['.tsv'])
    validate_file(phmmer, ['.tsv'])

    # Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

	# Set the path for the output file
    result_path = os.path.join(output_dir, f"{prefix}_hyp_table.csv")
    final_path = os.path.join(output_dir, f"{prefix}_hyp_final.csv")

    for path in [result_path, final_path]:
        if os.path.exists(path):
            response = input(f"Output file '{path}' already exists. Overwrite? (y/N): ").strip().lower()
            if response != 'y':
                print("Aborting to avoid overwriting the existing file.")
                sys.exit(1)

	#==========================================================#

	# Load and process blastp results
    blastp_columns = ['Query ID', 'Subject title']
    blastp_regex = r'UniRef90_\w+ (.+?) n=\d+ Tax=.+ RepID=\w+'
    blastp_df = load_and_extract_proteins(blastp, blastp_columns, 'Subject title', blastp_regex)
    blastp_results = collect_protein_results(blastp_df, 'Query ID', 'blastp_protein')

	# Load and process phmmer results
    phmmer_columns = ['target name', 'description of target']
    phmmer_regex = r'(.+?) n=\d+ Tax=.+ RepID=\w+'
    phmmer_df = load_and_extract_proteins(phmmer_columns, phmmer_columns, 'description of target', phmmer_regex)
    phmmer_results = collect_protein_results(phmmer_df, 'target name', 'phmmer_protein')

    # Merge results
    result_df = pd.merge(blastp_results, phmmer_results, left_index=True, right_index=True, how='outer')
    result_df = result_df.fillna(';'.join(['nan'] * MAX_PROTEINS))

    # Write results to files
    write_results(result_df, result_path, final_path)

if __name__ == '__main__':
	main()
