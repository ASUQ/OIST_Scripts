#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
	extract_tRNA -i <LOGFILE> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Purpose:
	Extract tRNA record from prokka log file.

Options:
    -h, --help                  Show this
    -i, --input		            Log file from prokka
    -p, --prefix <PREFIX>       Prefix for output file [default: tRNA]
    -o, --out <OUT>             Output directory path (string) [default: ./]

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
        description="Extract tRNA record from prokka log file",
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="prokka log file"
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
        default='tRNA',
		help="Prefix of output file [default: tRNA]"
	)

    return parser.parse_args()


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


def extract_tRNA_info(logfile_path):
    """
    Extract tRNA records from the prokka log file.

    Args:
        logfile_path (str): Path to the prokka log file.

    Returns:
        dict: Counts of each tRNA record.
    """

    # Dictionary to store the counts of each tRNA-XXX xxx combination
    tRNA_counts = {}

    # Flags to control reading lines between specific patterns
    in_target_section = False

    # Read the logfile
    with open(logfile_path, 'r') as logfile:
        for line in logfile:
            if line.startswith('[') and 'Running: aragorn -l' in line:
                in_target_section = True
            elif line.startswith('[') and 'Found' in line and 'tRNAs' in line:
                in_target_section = False
            elif in_target_section:
                if line.startswith('[') and 'tRNA-' in line:
                    # Extract XXX and xxx
                    parts = line.split()
                    tRNA_part = parts[2]
                    xxx_part = parts[-1].replace('t', 'u')

                    # Store the result in a tuple
                    key = tRNA_part + ' ' + xxx_part

                    # Increment the count for the extracted combination
                    if key in tRNA_counts:
                        tRNA_counts[key] += 1
                    else:
                        tRNA_counts[key] = 1

    return tRNA_counts

tRNA_list = [
	'tRNA-Phe (aaa)', 'tRNA-Phe (gaa)', 'tRNA-Leu (uaa)', 'tRNA-Leu (caa)', 'tRNA-Leu (aag)',
	'tRNA-Leu (gag)', 'tRNA-Leu (uag)', 'tRNA-Leu (cag)', 'tRNA-Ile (aau)', 'tRNA-Ile (gau)',
	'tRNA-Ile (uau)', 'tRNA-Met (cau)', 'tRNA-Val (aac)', 'tRNA-Val (gac)', 'tRNA-Val (uac)',
	'tRNA-Val (cac)', 'tRNA-Ser (aga)', 'tRNA-Ser (gga)', 'tRNA-Ser (uga)', 'tRNA-Ser (cga)',
	'tRNA-Pro (aag)', 'tRNA-Pro (ggg)', 'tRNA-Pro (ugg)', 'tRNA-Pro (cgg)', 'tRNA-Thr (agu)',
	'tRNA-Thr (ggu)', 'tRNA-Thr (ugu)', 'tRNA-Thr (cgu)', 'tRNA-Ala (agc)', 'tRNA-Ala (ggc)',
	'tRNA-Ala (ugc)', 'tRNA-Ala (cgc)', 'tRNA-Tyr (aua)', 'tRNA-Tyr (gua)', 'tRNA-Pyl (cua)',
	'tRNA-His (aug)', 'tRNA-His (gug)', 'tRNA-Gln (uug)', 'tRNA-Gln (cug)', 'tRNA-Asn (auu)',
	'tRNA-Asn (guu)', 'tRNA-Lys (uuu)', 'tRNA-Lys (cuu)', 'tRNA-Asp (auc)', 'tRNA-Asp (guc)',
	'tRNA-Glu (uuc)', 'tRNA-Glu (cuc)', 'tRNA-Cys (aca)', 'tRNA-Cys (gca)', 'tRNA-Trp (cca)',
	'tRNA-Arg (acg)', 'tRNA-Arg (gcg)', 'tRNA-Arg (ucg)', 'tRNA-Arg (ccg)', 'tRNA-Ser (acu)',
	'tRNA-Ser (gcu)', 'tRNA-Arg (ucu)', 'tRNA-Arg (ccu)', 'tRNA-Gly (acc)', 'tRNA-Gly (gcc)',
	'tRNA-Gly (ucc)', 'tRNA-Gly (ccc)']


def main():
    # Parse arguments
    args = parse_arguments()

	# Package import
    import os
    import sys

	# Extract argument values
    logfile_path = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.out)
    prefix = args.prefix

	# Check file existence and extension
    validate_file(logfile_path, ['.log'])

	# Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

	# Set the path for the output file
    output_file = os.path.join(output_dir, f"{prefix}.txt")
    if os.path.exists(output_file):
        response = input(f"Output file '{output_file}' already exists. Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("Aborting to avoid overwriting the existing file.")
            sys.exit(1)

	##===============================================================

	# Count tRNAs
    tRNA_counts = extract_tRNA_info(logfile_path)
    print(tRNA_counts)

    output = [(tRNA, tRNA_counts.get(tRNA, 0)) for tRNA in tRNA_list]

    # Return output
    with open(output_file, 'w') as output_handle:
        for i in range(len(output)):
            tRNA, num = output[i]
            output_handle.write(f'{tRNA},{num}')
            if i != len(output) - 1:
                output_handle.write('\n')


if __name__ == '__main__':
	main()
