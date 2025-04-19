#!/usr/bin/env python3
"""
Usage: (python) Busco_multigene_tree.py input_dir output_dir
"""
import sys
import os
import subprocess


# Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line


def run_mafft(input_file, output_file, cmd):
    """
    Runs the MAFFT alignment program on the specified file

    !!!Ensure MAFFT module is loaded on deigo before running!!!
    !!!Note that this is intended to be run on 64 threads. If this is not available, reduce the number!!!
    """
    print(f"Running MAFFT on {input_file}")
    mafftcmd = cmd.format(input_file, output_file)
    error_code = subprocess.check_call(mafftcmd, shell=True)

    if error_code != 0:
        print(
            f"Error occured while running mafft on {input_file}. Error code {error_code}"
        )


def run_trimal(input_file, output_file, cmd):
    """
    Runs the trimal alignment trimming on the specified file
    """
    print(f"Running trimal on {input_file}")
    trimcmd = cmd.format(input_file, output_file)
    error_code = subprocess.check_call(trimcmd, shell=True)

    if error_code != 0:
        print(
            f"Error occured while running trimal on {input_file}. Error code {error_code}"
        )


def run_iqtree(input_file, output_prefix, cmd):
    """
    Runs iqtree tree reconstruction on the specified file. Since the pipeline
    wont know the ideal models, it let's iqtree run modelfinder

    !!!Note that this is intended to be run on 64 threads. If this is not available, reduce the number!!!
    """
    print(f"Running iqtree on {input_file}")
    treecmd = cmd.format(input_file, output_prefix)
    error_code = subprocess.check_call(treecmd, shell=True)

    if error_code != 0:
        print(
            f"Error occured while running iqtree on {input_file}. Error code {error_code}"
        )


## Parse arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
fractions = sys.argv[3]
MAFFT_cmd = sys.argv[4]
trimal_cmd = sys.argv[5]
iqtree_cmd = sys.argv[6]

## Allow multiple fractions
if "," in fractions:
    fractions = list(map(float, fractions.split(",")))
else:
    fractions = [float(fractions)]


## Collect sequences by gene
gene_dict = {}  # gene_name: [list of orgs that have it]
org_list = []  # all organism names seen

for root, dirs, files in os.walk(input_dir):
    for file in files:
        # Look for single_copy_busco faa file
        if file.endswith(".faa") and root.endswith("single_copy_busco_sequences"):
            filepath = os.path.join(root, file)
            org_name = root.replace(input_dir, "").strip("/").split("/")[0]
            gene_name = file.strip(".faa")

            if org_name not in org_list:
                org_list.append(org_name)

            if gene_name in gene_dict:
                gene_dict[gene_name].append(org_name)
            else:
                gene_dict[gene_name] = [org_name]

            # Create output file
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)

            # Write the busco sequences to their respective folders
            with open(os.path.join(output_dir, gene_name + ".faa"), "a") as out:
                for line in gen_line_reader(os.path.join(root, file)):
                    # Change header name to organism name
                    if line.strip():
                        if line.startswith(">"):
                            out.write(">" + org_name + "\n")
                        else:
                            out.write(line)


## Selecting "shared" genes by fraction
fract_dict = {}  # fraction: [genes meeting the threshold]
for fraction in fractions:
    shared_genes = []
    threshold = len(org_list) * fraction
    for gene in gene_dict:
        if len(gene_dict[gene]) >= threshold:
            shared_genes.append(gene)
        fract_dict[fraction] = shared_genes

# write out which genes were picked at each fraction
for fraction in fractions:
    with open(
        os.path.join(output_dir, f"fraction{fraction}_analyzed_genes.txt"), "w"
    ) as gene_out:
        gene_out.write(
            f"Number of genes considered: {len(fract_dict[fraction])}\nAnalyzed genes:\n"
        )
        for g in fract_dict[fraction]:
            gene_out.write(g + "\n")


## Align and Trim the "smallest" fraction
smallest = min(fractions)

for gene in fract_dict[smallest]:
    run_mafft(
        os.path.join(output_dir, gene + ".faa"),
        os.path.join(output_dir, gene + ".faa_aligned"),
        MAFFT_cmd,
    )
    run_trimal(
        os.path.join(output_dir, gene + ".faa_aligned"),
        os.path.join(output_dir, gene + ".faa_aligned_trimmed"),
        trimal_cmd,
    )


## Concatenate alignments and Build Trees
for fraction in fractions:
    concat_dict = {}
    for gene in fract_dict[fraction]:
        for line in gen_line_reader(
            os.path.join(output_dir, gene + ".faa_aligned_trimmed")
        ):
            if line.startswith(">"):
                org_name = line.strip(">").strip()
                if org_name not in concat_dict:
                    concat_dict[org_name] = {}
                concat_dict[org_name][gene] = []
            else:
                concat_dict[org_name][gene].append(line.strip())

    # Prepare dictionary for outputting
    # first concatenate the lists of sequences into a single string
    for org in concat_dict:
        for gene in concat_dict[org]:
            concat_dict[org][gene] = "".join(concat_dict[org][gene])
    # Second, find the genes that are missing, and fill them with dashes, the same lenght as the strings
    # This is to allow iqtree to use the concatenated sequences, since they need to be the same length
    for org in concat_dict:
        for gene in fract_dict[fraction]:
            if gene not in concat_dict[org]:
                # find the required length of dashes and add them to the dictionary
                for org2 in concat_dict:
                    if gene in concat_dict[org2]:
                        concat_dict[org][gene] = "-" * len(concat_dict[org2][gene])

    # Write concatenated file
    with open(
        os.path.join(
            output_dir, f"fraction{fraction}_concatenated_trimmed_alignment.faa"
        ),
        "w",
    ) as concat_file:
        for org in concat_dict:
            concat_file.write(">" + org + "\n")
            for gene in fract_dict[
                fraction
            ]:  # Use shared genes here to ensure the order remains the same
                concat_file.write("".join(concat_dict[org][gene]))
            concat_file.write("\n")

    # Lastly, run iqtree on the concatenated file (per given fraction)
    run_iqtree(
        os.path.join(
            output_dir, f"fraction{fraction}_concatenated_trimmed_alignment.faa"
        ),
        os.path.join(output_dir, f"fraction{fraction}_iqtree_output"),
        iqtree_cmd,
    )
