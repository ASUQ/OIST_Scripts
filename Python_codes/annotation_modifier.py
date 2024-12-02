#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    annotation_modifier -g <GFF> -s <SCAFFOLDS> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Purpose:
	Modify the locus name of gff file to make it compatible to Tablet

Options:
    -h, --help                      Show this
    -g, --gff <GFF>                 GFF file of annotaion (string) - generated from prokka
    -s, --scaffolds <SCAFFOLDS>     Fasta file of annotated scaffolds
    -p, --prefix <PREFIX>           Prefix for output file. [default: modified_<gfffile>]
    -o, --out <OUT>                 Output directory path (string) [default: ./]

Required Packages:
    Biopython, gffutils

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
        description= "Modify the locus name of gff file to make it compatible to Tablet",
        epilog="Required package: Biopython, gffutils",
        formatter_class=argparse.RawDescriptionHelpFormatter
	)

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        "-g", "--gff",
        type=str,
        required=True,
        help="GFF file of annotaion (string) - generated from prokka"
	)

    required.add_argument(
        "-s", "--scaffolds",
		type=str,
		required=True,
		help="Fasta file of annotated scaffolds"
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
		help="Prefix of output file [default: modified_<gfffile>]"
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


def extract_scaffold_ids(scaffolds_file):
    """
    Extract scaffold IDs from a FASTA file.

    Args:
        scaffolds_file (str): Path to the FASTA file.

    Returns:
        list: List of scaffold IDs.
    """
    from Bio import SeqIO

    scaffolds_ID = []
    with open(scaffolds_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            scaffolds_ID.append(record.id)
    return scaffolds_ID


def replace_locus_ids(gff_file, output_tmp, locus_to_scaffolds):
    """
    Replace LOCUS IDs in the GFF file with scaffold IDs.

    Args:
        gff_file (str): Path to the input GFF file.
        output_tmp (str): Path to the temporary output file.
        locus_to_scaffolds (dict): Mapping of LOCUS IDs to scaffold IDs.
    """
    comments = []
    with open(gff_file, 'r') as gff, open(output_tmp, 'w') as tmp:
        while True:
            line = gff.readline()

            for key, value in locus_to_scaffolds.items():
                line = line.replace(key, value)
            tmp.write(line)

            if line == '##FASTA\n':
                break
            elif line[:2] == '##':
                comments.append(line)
    return comments


def modify_features(tmp_file, output_file):
    """
    Modify features in the GFF file.

    Args:
        tmp_file (str): Path to the temporary GFF file.
        output_file (str): Path to the final GFF file.
    """
    import gffutils

	# Create a GFF database
    db = gffutils.create_db(tmp_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy="merge")

	# Create a list to store all features with modifications
    all_features = []

    for feature in db.all_features():
        if feature.featuretype in ['CDS', 'tRNA', 'rRNA'] and 'Name' not in feature.attributes and 'product' in feature.attributes:
            # Create a modified attribute dictionary
            updated_attributes = dict(feature.attributes)
            updated_attributes['Name'] = feature.attributes['product']

			# Format the modified feature as a GFF string
            attr_str = ';'.join([f"{key}={value[0]}" for key, value in updated_attributes.items()])
            modified_feature_str = '\t'.join([
                feature.seqid,
                feature.source,
                feature.featuretype,
                str(feature.start),
                str(feature.end),
                str(feature.score),
                feature.strand,
                feature.frame,
                attr_str
            ])
            all_features.append(modified_feature_str)
        else:
            # For unmodified features or those not meeting modification criteria, add their existing GFF strings
            all_features.append(str(feature))

	# Write all features to the updated GFF file
    with open(output_file, 'a') as modified_gff:
        for feature_str in all_features:
            modified_gff.write(feature_str + '\n')


def append_fasta(output_file, scaffolds_file):
    """
    Append the FASTA section to the GFF file.

    Args:
        output_file (str): Path to the final GFF file.
        scaffolds_file (str): Path to the FASTA file.
    """
    with open(output_file, 'a') as modified_gff:
        modified_gff.write('##FASTA\n')
        with open(scaffolds_file, 'r') as fasta:
            for line in fasta:
                modified_gff.write(line)


def main():
    # Parse arguments
    args = parse_arguments()

    # Check the installed packages
    required_packages = ['Bio', 'gffutils']
    check_required_packages(required_packages)

	# Package import
    from Bio import SeqIO
    import gffutils
    import os
    import sys

    # Extract argument values
    gff_file = args.gff
    scaffolds = args.scaffolds
    output_dir = args.out
    prefix = args.prefix

	# Check file existence and extension
    validate_file(gff_file, ['.gff'])
    validate_file(scaffolds, ['.fasta', '.fas', '.fa', '.fna'])

	# Make new directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Check prefix name
    if prefix is None:
        prefix = 'modified_' + str(os.path.splitext(os.path.basename(gff_file))[0])

	# Set the path for the output file
    output_file = os.path.join(output_dir, f"{prefix}.gff")
    if os.path.exists(output_file):
        response = input(f"Output file '{output_file}' already exists. Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("Aborting to avoid overwriting the existing file.")
            sys.exit(1)

	# Temporary file
    output_tmp = os.path.join(output_dir, f"{prefix}.gff.tmp")

    ##======================================================================================

	## Step 1: Replace LOCUS IDs
    scaffolds_ID = extract_scaffold_ids(scaffolds)
    locus_to_scaffolds = {}
    for i in range(len(scaffolds_ID)):
        locus_to_scaffolds[f'gnl|Prokka|LOCUS_{i+1}'] = scaffolds_ID[i]

    # To avoid mismatch (e.g. LOCUS_1 match with LOCUS_11)
    locus_to_scaffolds = dict(reversed(locus_to_scaffolds.items()))

    comments = replace_locus_ids(gff_file, output_tmp, locus_to_scaffolds)

    # Store comments into modified file
    with open(output_file, 'w') as modified:
        modified.writelines(comments)


	## Step 2: Modify features
    modify_features(output_tmp, output_file)


    ## Step 3: Append Fasta
    append_fasta(output_file, scaffolds)


	## Step 4: Cleanup temporary file
    os.remove(output_file + '.tmp')


if __name__ == '__main__':
	main()
