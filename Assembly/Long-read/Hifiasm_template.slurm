#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=200G
#SBATCH --cpus-per-task 32
#SBATCH --nodes 1
#SBATCH --job-name=MMDD_Job
#SBATCH --output=/path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.name@oist.jp

### Slurm information: https://groups.oist.jp/scs/use-slurm
### Save as /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm

#===============================================================================================#

# Description:
# This script is for assembling pacbio hifi data with Deigo (computer cluster at OIST)
#
# Workflow:
# 1. Assemble long read with Hifiasm
# 2. Annotate scaffolds taxonomically with megablast
# 3. Map reads onto assembly with minimap2 and samtools
# 4. Summarize results with Blobtools
#
# Expected Inputs:
# - Pacbio hifi read
#
# Required Softwares:
# - Hifiasm
# - Standalone BLAST
# - minimap2
# - Samtools
# - Blobtools v1
# Note: Miniconda needs to be installed.
#
# Usage:
# - The paths of inputs and softwares, and parameters are hardcoded.
# - Make sure to update the script with the correct paths and input filenames before execution.

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
HIFIASM='<Hifiasm directory>'
MINIMAP2='<minimap2 directory>'
SAMTOOLS='<samtools directory>'

# Database paths
BLOBTOOLS_DB='Path to blobtools nodesDB'
BLAST_NT_DB='Path to nt database from ncbi'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'

SAMPLE='<sample name>'

long_read='<Long read>.fastq.gz'	# Adjust the minimap2 option based on your read types
#------------------------------------#

#-------------Parameters-------------#
THREADS='<number of CPUs>'	# More than 16 is recommended
# MEMORY='<size of RAM>'
#------------------------------------#


#===============================================================================================#

# DON'T MODIFY THE SCRIPT BELOW!!

set_software_paths() {
	# blast
	BLASTN="${BLAST}/bin/blastn"

	# blobtools
	BLOBTOOLS="${BLOBTOOLS}/blobtools"

	# hifiasm
	HIFIASM="${HIFIASM}/hifiasm"

	# minimap2
	MINIMAP2="${MINIMAP2}/minimap2"

	# samtools
	SAMTOOLS="${SAMTOOLS}/samtools"

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

run_hifiasm() {
    echo '-----------------------------------------------'
	echo 'Running Hifiasm'
	echo "$(date)"
	echo '-----------------------------------------------'

	$HIFIASM -o "${RUN_DIR}/${SAMPLE}" -t $((THREADS)) "${long_read}"

    if [ $? -ne 0 ]; then
        echo "Hifiasm failed" >&2
        exit 1
    fi

	# Convert gfa file into fasta file
	awk '/^S/{print ">"$2;print $3}' "${RUN_DIR}/${SAMPLE}.bp.p_ctg.gfa" > "${RUN_DIR}/assembly.fasta"

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

    ${BLASTN} -task megablast -query "${RUN_DIR}/assembly.fasta" \
		-db "${BLAST_NT_DB}" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
		-culling_limit 5 -num_threads $((THREADS)) -evalue 1e-25 -max_target_seqs 5 > "${RUN_DIR}/assembly_vs_nt.blast"

    if [ $? -ne 0 ]; then
        echo "megablast failed" >&2
        exit 1
    fi

	return 0
}

run_minimap2() {
	echo '-----------------------------------------------'
	echo 'Running minimap2'
	echo "$(date)"
	echo '-----------------------------------------------'

	mkdir -p minimap2

	${MINIMAP2} -a -x map-hifi -t $((THREADS)) "${RUN_DIR}/assembly.fasta" "${long_read}" > "${RUN_DIR}/minimap2/${SAMPLE}.sam"

    if [ $? -ne 0 ]; then
        echo "minimap2 failed" >&2
        exit 1
    fi

	return 0
}

run_samtools() {
	echo '-----------------------------------------------'
	echo 'Running samtools'
	echo "$(date)"
	echo '-----------------------------------------------'

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}.sam" -o "${RUN_DIR}/minimap2/${SAMPLE}.sorted.bam"

	rm -f "${RUN_DIR}/minimap2/${SAMPLE}.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}.sorted.bam"

	if [ $? -ne 0 ]; then
        echo "samtools failed" >&2
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

    ${BLOBTOOLS} create -i "${RUN_DIR}/assembly.fasta" -b "${RUN_DIR}/minimap2/${SAMPLE}.sorted.bam" \
		-t "${RUN_DIR}/assembly_vs_nt.blast" -o "${RUN_DIR}/blobtools" --db "${BLOBTOOLS_DB}"

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

	## Run Hifiasm to assemble the reads
    run_hifiasm

	## Run megablast to taxonomically annotate the scaffolds
    run_blastn

	## Run minimap2 to map the read to the reference
	run_minimap2

	## Run samtools to convert sam file to sorted bam file and index
	run_samtools

	## Run blobtools to summarize the blastn result
    run_blobtools

    echo 'JOB COMPLETED'

	return 0
}


# Execute main function
main

### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
