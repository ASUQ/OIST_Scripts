#!/home/a/akito-shima/miniconda3/envs/core/bin/python
# -*- coding: utf-8 -*-

"""
Usage:
    sequence_extractor -i <INPUT> -f <FASTA> [-p <PREFIX>] [-o <OUT>] [-h|--help]

Options:
    -h, --help                  Show this
    -i, --input                 List of seqeunces you want (string) - list by "\n"
    -f, --fasta <FASTA>         Fasta file of sequences (string) - generated from genome assembler
    -p, --prefix <PREFIX>       Prefix for output file.
    -o, --out <OUT>             Output directory path (string) [default: ./]

Required Packages:
    Biopython

Author: Akito Shima (ASUQ)
Email: akito-shima@oist.jp
"""

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    import os.path
    import sys


    ## Command input using argparse
    parser = argparse.ArgumentParser(
                                    epilog="required package: Biopython"
                                    )

    parser.add_argument("-i", "--input",
                        type=str,
                        required=True,
                        help="List of sequences you want (string) - list by LF"
                        )

    parser.add_argument("-f", "--fasta",
                        type=str,
                        required=True,
                        help="Fasta file of sequences (string) - generated from genome assembler"
                        )

    parser.add_argument("-o", "--out",
                        type=str,
                        default='./',
                        help="Output directory path  [default: ./]"
                        )

    parser.add_argument("-p", "--prefix",
                        type=str,
                        help="Output directory path  [default: input file name]"
                        )

    # Parse the inputs from command-line
    args = parser.parse_args()
    nodes = args.input
    sequences = args.fasta
    out = args.out
    prefix = args.prefix

    # Check file existence and extension
    with open(nodes, 'r') as f:
        root, ext = os.path.splitext(nodes)
        if ext.lower() != '.txt':
            raise ValueError(f"File '{nodes}' has an incorrect extension. Expected '.txt' file.")

        # Store the name of requested sequences in list
        names = []
        for line in f:
            if len(line.strip()) > 0:
                names.append(line.strip())

        # Remove duplicates and store as immutable tuple
        names = list(dict.fromkeys(names))
        names = tuple(names)

    # Check sequences file existence and extension
    with open(sequences, 'r') as f:
        root, ext = os.path.splitext(sequences)
        if ext.lower() not in ['.fasta', '.fas', '.fa', '.fna']:
            raise ValueError(f"File '{sequences}' has an incorrect extension. Expected fasta file.")


    # Check output directory path
    if out[-1] != '/':
        out += '/'
    # Make new directory if directory does not exist
    if not os.path.isdir(out):
        os.makedirs(out)

    # Check prefix name
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(nodes))[0]

    # Extract sequences of the targets and save as fasta file
    sequence_records = dict()
    with open(sequences, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in names:
                sequence_records[record.id] = record

    # print(sequence_records)

    # print(f'sequence_records:{sequence_records}')
    if len(sequence_records):
        output = []
        for name in names:
            output.append(sequence_records[name])
        SeqIO.write(output, out + prefix + '.fasta', 'fasta')
    else:
        print('No sequence was found.')


    # # Terminate the process
    # sys.exit()
