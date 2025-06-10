#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=500G
#SBATCH --cpus-per-task 128
#SBATCH --nodes 1
#SBATCH --job-name=MMDD_Job
#SBATCH --output=/path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.name@oist.jp

### Slurm information: https://groups.oist.jp/scs/use-slurm
### Save as /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm

#===============================================================================================#

# Description:
# This script is for assembling Illumina short-read data with Deigo (computer cluster at OIST).
#
# Workflow:
# 1. Assemble short-reads with SPAdes
# 2. Annotate scaffolds taxonomically with megablast
# 3. Summarize results with Blobtools
#
# Expected Inputs:
# - Illumina pair-end reads (quality controlled)
#
# Required Softwares:
# - SPAdes
# - Standalone BLAST
# - Blobtools v1
# Note: Miniconda needs to be installed.
#
# Usage:
# - The paths of inputs and softwares, and parameters are hardcoded.
# - Make sure to update the script with the correct paths and input filenames before execution.
# - For metagenome assembly, add --meta option and remove --careful option in SPAdes

#===============================================================================================#

set -euo pipefail

#===============================================================================================#

## Edit this part before running the script

#-------------Softwares--------------#
# Conda env paths
MINICONDA_DIR='Path to miniconda directory'
BLOBTOOLS_CONDA='Path to blobtools conda environment'

# Software paths
BLAST='<BLAST+ directory>'
BLOBTOOLS='<Blobtools directory>'
SPADES='<SPAdes directory>'		# Version 4 is recommended

# Database paths
BLOBTOOLS_DB='Path to blobtools nodesDB'
BLAST_NT_DB='Path to nt database from ncbi'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'

read_f='<Illumina forward read>.fastq.gz'
read_r='<Illumina reverse read>.fastq.gz'
#------------------------------------#

#-------------Parameters-------------#
THREADS="$SLURM_CPUS_PER_TASK"
MEMORY="$SLURM_MEM_PER_NODE"
#------------------------------------#


#===============================================================================================#

# DON'T MODIFY THE SCRIPT BELOW!!

set_software_paths() {
	# blast
	BLASTN="${BLAST}/bin/blastn"

	# blobtools
	BLOBTOOLS="${BLOBTOOLS}/blobtools"

	# spades
	SPADES="${SPADES}/bin/spades.py"

	return 0
}

activate_conda() {
    source "${MINICONDA_DIR}/etc/profile.d/conda.sh"
    if ! command -v conda &> /dev/null; then
        echo "conda could not be found"
        exit 1
    fi

	return 0
}

run_spades() {
    echo '-----------------------------------------------'
	echo 'Running SPAdes'
	echo "$(date)"
	echo '-----------------------------------------------'

	# module installation
	module purge
	module load python/3.11.4	# Default python version on Deigo is v3.6.8, v3.8+ is needed for SPAdes v4.0+

	# For metagenome assembly, add --meta option and remove --careful option in SPAdes
    ${SPADES} -o "${RUN_DIR}" -1 "${read_f}" -2 "${read_r}" \
		-k 21,33,55,77,99,119,127 --careful --threads $((THREADS)) --memory $((MEMORY))

    if [ $? -ne 0 ]; then
        echo "SPAdes failed" >&2
        exit 1
    fi

	module purge

	return 0
}

run_blastn() {
	echo '-----------------------------------------------'
	echo 'Running megablast'
	echo "$(date)"
	echo '-----------------------------------------------'

	# To avoid warnings
	BLASTDB="$BLAST_NT_DB/.."
	export BLASTDB

    ${BLASTN} -task megablast -query "${RUN_DIR}/scaffolds.fasta" \
		-db "${BLAST_NT_DB}" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
		-culling_limit 5 -num_threads $((THREADS)) -evalue 1e-25 -max_target_seqs 5 > "${RUN_DIR}/scaffolds_vs_nt.blast"

    if [ $? -ne 0 ]; then
        echo "megablast failed" >&2
        exit 1
    fi

	return 0
}

run_blobtools() {
	echo '-----------------------------------------------'
	echo 'Running blobtools'
	echo "$(date)"
	echo '-----------------------------------------------'

    conda activate "${BLOBTOOLS_CONDA}"

	if [ $? -ne 0 ]; then
        echo "Failed to activate blobtools conda environment" >&2
        exit 1
    fi

    ${BLOBTOOLS} create -i "${RUN_DIR}/scaffolds.fasta" -y spades -t "${RUN_DIR}/scaffolds_vs_nt.blast" \
		-o "${RUN_DIR}/blobtools" --db "${BLOBTOOLS_DB}"

    ${BLOBTOOLS} plot -i "${RUN_DIR}/blobtools.blobDB.json" -r superkingdom -o "${RUN_DIR}/"
    ${BLOBTOOLS} plot -i "${RUN_DIR}/blobtools.blobDB.json" -r phylum -o "${RUN_DIR}/"
    ${BLOBTOOLS} plot -i "${RUN_DIR}/blobtools.blobDB.json" -r order -o "${RUN_DIR}/"


    ${BLOBTOOLS} view -i "${RUN_DIR}/blobtools.blobDB.json" -r all -o "${RUN_DIR}/"

	# Remove comments and extract data from blobtools view output
	echo 'name,length,GC,N,coverage,superkingdom,phylum,order,family,genus,species' > "${RUN_DIR}/blobtools.blobDB.table.csv"
	grep -v '^##' "${RUN_DIR}/blobtools.blobDB.table.txt" | \
		grep -v '^#' |
		cut --complement -f7,8,10,11,13,14,16,17,19,20,22,23 |
		tr -s '\t' ',' >> "${RUN_DIR}/blobtools.blobDB.table.csv"

	conda deactivate

	return 0
}

main() {
    ## Add conda to your PATH
	activate_conda

	## Create running directory if doesn't exist
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}" || exit

	## Set software paths
	set_software_paths

	## Run SPAdes to assemble the reads
    run_spades

	## Run megablast to taxonomically annotate the scaffolds
    run_blastn

	## Run blobtools to summarize the blastn result
    run_blobtools

    echo 'JOB COMPLETED'

	return 0
}


# Execute main function
main

### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
