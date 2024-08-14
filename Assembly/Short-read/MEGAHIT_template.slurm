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
# This script is for assembling Illumina short-read data with Deigo (computer cluster at OIST).
#
# Workflow:
# 1. Assemble short-reads with MEGAHIT
# 2. Annotate scaffolds taxonomically with megablast
# 3. Map reads onto assembly with bowtie2 and samtools
# 4. Summarize results with Blobtools
#
# Expected Inputs:
# - Illumina pair-end reads (quality controlled)
#
# Required Softwares:
# - MEGAHIT
# - Standalone BLAST
# - bowtie2
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
MEGAHIT_CONDA='Path to megahit conda environment'

# Software paths
BLAST='<BLAST+ directory>'
BLOBTOOLS='<Blobtools directory>'
BOWTIE2='<bowtie2 directory>'
SAMTOOLS='<samtools directory>'

# Database paths
BLOBTOOLS_DB='Path to blobtools nodesDB'
BLAST_NT_DB='Path to nt database from ncbi'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'

SAMPLE='<sample name>'

read_f='<Illumina forward read>.fastq.gz'
read_r='<Illumina reverse read>.fastq.gz'
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

	# bowtie
	BOWTIE2_INDEX="${BOWTIE2}/bowtie2-build"
	BOWTIE2="${BOWTIE2}/bowtie2"

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

run_megahit() {
    echo '-----------------------------------------------'
	echo 'Running MEGAHIT'
	echo "$(date)"
	echo '-----------------------------------------------'

	conda activate "${MEGAHIT_CONDA}"

	if [ $? -ne 0 ]; then
        echo "Failed to activate megahit conda environment" >&2
        exit 1
    fi

	megahit -1 "${read_f}" -2 "${read_r}" -o "${RUN_DIR}/Megahit" --out-prefix "${SAMPLE}_Megahit" \
		--no-mercy --num-cpu-threads $((THREADS)) --memory 1.0

    if [ $? -ne 0 ]; then
        echo "MEGAHIT failed" >&2
        exit 1
    fi

	conda deactivate

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

    ${BLASTN} -task megablast -query "${RUN_DIR}/Megahit/${SAMPLE}_Megahit.contigs.fa" \
		-db "${BLAST_NT_DB}" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
		-culling_limit 5 -num_threads $((THREADS)) -evalue 1e-25 -max_target_seqs 5 > "${RUN_DIR}/Megahit/contigs_vs_nt.blast"

    if [ $? -ne 0 ]; then
        echo "megablast failed" >&2
        exit 1
    fi

	return 0
}

run_bowtie2() {
	echo '-----------------------------------------------'
	echo 'Running bowtie2'
	echo "$(date)"
	echo '-----------------------------------------------'

	mkdir -p Megahit/bowtie2

	${BOWTIE2_INDEX} -f -q --threads $((THREADS)) "${RUN_DIR}/Megahit/${SAMPLE}_Megahit.contigs.fa" "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit"
	${BOWTIE2} -q --reorder --threads $((THREADS)) --time --met-stderr --met 10 \
		-x "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit" -1 "${read_f}" -2 "${read_r}" > "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit".sam

    if [ $? -ne 0 ]; then
        echo "bowtie2 failed" >&2
        exit 1
    fi

	rm -f "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.*.bt2"

	return 0
}

run_samtools() {
	echo '-----------------------------------------------'
	echo 'Running samtools'
	echo "$(date)"
	echo '-----------------------------------------------'

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.sam" -o "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.sorted.bam"

	rm -f "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.sorted.bam"

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

    ${BLOBTOOLS} create -i "${RUN_DIR}/Megahit/${SAMPLE}_Megahit.contigs.fa" -b "${RUN_DIR}/Megahit/bowtie2/${SAMPLE}_Megahit.sorted.bam" \
		-t "${RUN_DIR}/Megahit/contigs_vs_nt.blast" -o "${RUN_DIR}/Megahit/blobtools" --db "${BLOBTOOLS_DB}"

    ${BLOBTOOLS} plot -i "${RUN_DIR}/Megahit/blobtools.blobDB.json" -r superkingdom -o "${RUN_DIR}/Megahit/"
    ${BLOBTOOLS} plot -i "${RUN_DIR}/Megahit/blobtools.blobDB.json" -r phylum -o "${RUN_DIR}/Megahit/"
    ${BLOBTOOLS} plot -i "${RUN_DIR}/Megahit/blobtools.blobDB.json" -r order -o "${RUN_DIR}/Megahit/"


    ${BLOBTOOLS} view -i "${RUN_DIR}/Megahit/blobtools.blobDB.json" -r all -o "${RUN_DIR}/Megahit/"

	# Remove comments and extract data from blobtools view output
	echo 'name,length,GC,N,coverage,superkingdom,phylum,order,family,genus,species' > "${RUN_DIR}/Megahit/blobtools.blobDB.table.csv"
	grep -v '^##' "${RUN_DIR}/Megahit/blobtools.blobDB.table.txt" | \
		grep -v '^#' |
		cut --complement -f7,8,10,11,13,14,16,17,19,20,22,23 |
		tr -s '\t' ',' >> "${RUN_DIR}/Megahit/blobtools.blobDB.table.csv"

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

	## Run MEGAHIT to assemble the reads
    run_megahit

	## Run megablast to taxonomically annotate the scaffolds
    run_blastn

	## Run bowtie2 to map the read to the reference
	run_bowtie2

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
