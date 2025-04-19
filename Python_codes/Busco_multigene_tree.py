#!/usr/bin/env python3
"""
Usage: (python) Busco_multigene_tree.py input_dir output_dir
"""
import sys
import os
import subprocess

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

def run_mafft (input_file, output_file, cmd):
    """
    Runs the MAFFT alignment program on the specified file
    
    !!!Ensure MAFFT module is loaded on deigo before running!!!
    !!!Note that this is intended to be run on 64 threads. If this is not available, reduce the number!!!
    """
    print("Running MAFFT on {}".format(input_file))
    mafftcmd = cmd.format(input_file,output_file)
    error_code = subprocess.check_call(mafftcmd, shell = True)
    #check if program ran correctly
    if error_code != 0:
        print("Error occured while running mafft on {}. Error code {}"\
        .format(input_file, error_code))
        
def run_trimal (input_file, output_file, cmd):
    """
    Runs the trimal alignment trimming on the specified file
    """
    print("Running trimal on {}".format(input_file))
    trimcmd = cmd.format(input_file,output_file)
    error_code = subprocess.check_call(trimcmd, shell = True)
    #check if program ran correctly
    if error_code != 0:
        print("Error occured while running trimal on {}. Error code {}"\
        .format(input_file, error_code))
        
def run_iqtree (input_file, output_prefix, cmd):
    """
    Runs iqtree tree reconstruction on the specified file. Since the pipeline
    wont know the ideal models, it let's iqtree run modelfinder
    
    !!!Note that this is intended to be run on 64 threads. If this is not available, reduce the number!!!
    """
    print("Running iqtree on {}".format(input_file))
    treecmd = cmd.format(input_file, output_prefix)
    error_code = subprocess.check_call(treecmd, shell = True)
    #check if program ran correctly
    if error_code != 0:
        print("Error occured while running iqtree on {}. Error code {}"\
        .format(input_file, error_code))
            
input_dir = sys.argv[1]
output_dir = sys.argv[2]
fractions = sys.argv[3]
MAFFT_cmd = sys.argv[4]
trimal_cmd = sys.argv[5]
iqtree_cmd = sys.argv[6]
if "," in fractions:
    fractions = fractions.split(",")
else:
    fractions = [fractions]
gene_dict = {}
org_list = []
for root, dirs, files in os.walk(input_dir): #iterate over all directories and subdirectories
    for file in files: #iterate over all files, only use the .faa files generated for single copy busco sequences
        if file.endswith(".faa"):
            if root.endswith("single_copy_busco_sequences"):
                #Parse important values
                filepath = os.path.join(root, file)
                org_name = root.replace(input_dir, "").strip("/").split("/")[0]
                if org_name not in org_list:
                    org_list.append(org_name)
                gene_name = file.strip(".faa")
                if gene_name in gene_dict:
                    gene_dict[gene_name].append(org_name)
                else:
                    gene_dict[gene_name] = [org_name]
                #Check if output folder exists, if not, create items
                if not os.path.isdir(output_dir):
                    os.mkdir(output_dir)
                #Write the busco sequences to their respective folders
                out = open(os.path.join(output_dir, gene_name + ".faa"),"a")
                for line in gen_line_reader(os.path.join(root, file)):
                    #Change header name to organism name
                    if line.strip():
                        if line.startswith(">"):
                            out.write(">" + org_name + "\n")
                        else:
                            out.write(line)
                out.close()

#Check which genes are shared by the required fraction of organisms
fract_dict = {}
for fraction in fractions:
    shared_genes = []
    for gene in gene_dict:
        if len(gene_dict[gene]) >= float(fraction) * len(org_list):
            shared_genes.append(gene)
        fract_dict[fraction] = shared_genes

#For reference, output the used genes
for fraction in fractions:
    gene_out = open(os.path.join(output_dir, "fraction{}_analyzed_genes.txt".format(fraction)), "w")
    gene_out.write("Number of genes considered: {}\nAnalyzed genes:\n".format(len(fract_dict[fraction])))
    for gene in fract_dict[fraction]:
        gene_out.write(gene +"\n")
    gene_out.close()

#Find the lowest fraction in the list (the one with the most genes), and run the alignment and trimming for the genes in this fraction
smallest = 100
for fraction in fractions:
    if float(fraction) < float(smallest):
        smallest = fraction
    

for gene in fract_dict[smallest]:
    #run the tools for alignment and trimming
    run_mafft(os.path.join(output_dir, gene + ".faa"),os.path.join(output_dir, gene + ".faa_aligned"), MAFFT_cmd)
    run_trimal(os.path.join(output_dir, gene + ".faa_aligned"), os.path.join(output_dir, gene + ".faa_aligned_trimmed"), trimal_cmd)
#Output the concatenated sequences per fraction in multiple steps
#Collect sequences for concatenating
for fraction in fractions:
    concat_dict = {}
    for gene in fract_dict[fraction]:
        for line in gen_line_reader(os.path.join(output_dir, gene + ".faa_aligned_trimmed")):
            if line.startswith(">"):
                org_name = line.strip(">").strip()
                if org_name not in concat_dict:
                    concat_dict[org_name] = {}
                concat_dict[org_name][gene] = []
            else:
                concat_dict[org_name][gene].append(line.strip())

    #Prepare dictionary for outputting
    #first concatenate the lists of sequences into a single string
    for org in concat_dict:
        for gene in concat_dict[org]:
            concat_dict[org][gene] = "".join(concat_dict[org][gene])
    #Second, find the genes that are missing, and fill them with dashes, the same lenght as the strings
    #This is to allow iqtree to use the concatenated sequences, since they need to be the same length
    for org in concat_dict:
        for gene in fract_dict[fraction]:
            if gene not in concat_dict[org]:
                #find the required length of dashes and add them to the dictionary
                for org2 in concat_dict:
                    if gene in concat_dict[org2]:
                       concat_dict[org][gene] = "-" * len(concat_dict[org2][gene]) 

    #Write concatenated file    
    concat_file = open(os.path.join(output_dir, "fraction{}_concatenated_trimmed_alignment.faa".format(fraction)), "w")    
    for org in concat_dict:
        concat_file.write(">" + org + "\n")
        for gene in fract_dict[fraction]: #Use shared genes here to ensure the order remains the same
            concat_file.write("".join(concat_dict[org][gene]))
        concat_file.write("\n")    
    concat_file.close()

    #Lastly, run iqtree on the concatenated file (per given fraction)
    run_iqtree(os.path.join(output_dir, "fraction{}_concatenated_trimmed_alignment.faa".format(fraction)), os.path.join(output_dir, "fraction{}_iqtree_output".format(fraction)), iqtree_cmd)    

              