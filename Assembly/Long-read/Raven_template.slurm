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
# This script is for assembling long-read data with Deigo (computer cluster at OIST)
#
# Workflow:
# 1. Assemble long read with Raven
# 2. Annotate scaffolds taxonomically with megablast
# 3. Map reads onto assembly with minimap2 and samtools
# 4. Summarize results with Blobtools
#
# Expected Inputs:
# - Pacbio or Oxford Nanopore long read    # You need to adjust the Flye option based on the read type
#
# Required Softwares:
# - Raven
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

set -euxo pipefail

#===============================================================================================#

## Edit this part before running the script

#-------------Softwares--------------#
# Conda env paths
MINICONDA_DIR='Path to miniconda directory'
BLOBTOOLS_CONDA='Path to blobtools conda environment'
RAVEN_CONDA='Path to raven conda environment'

# Software paths
BLAST='<BLAST+ directory>'
BLOBTOOLS='<Blobtools directory>'
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

run_raven() {
    echo '-----------------------------------------------'
	echo 'Running Raven'
	echo '-----------------------------------------------'

	conda activate "${RAVEN_CONDA}"

	if [ $? -ne 0 ]; then
        echo "Failed to activate raven conda environment" >&2
        exit 1
    fi

	raven --threads $((THREADS)) --graphical-fragment-assembly "${RUN_DIR}/assembly.gfa" "${long_read}" > "${RUN_DIR}/assembly.fasta"

    if [ $? -ne 0 ]; then
        echo "Raven failed" >&2
        exit 1
    fi

	conda deactivate

	return 0
}

run_blastn() {
	echo '-----------------------------------------------'
	echo 'Running megablast'
	echo '-----------------------------------------------'

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
	echo '-----------------------------------------------'

	mkdir -p minimap2

	# Adjust the read option based on your data
	# - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
	# - map-hifi - PacBio HiFi reads vs reference mapping

	${MINIMAP2} -a -x map-ont -t $((THREADS)) "${RUN_DIR}/assembly.fasta" "${long_read}" > "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sam"

    if [ $? -ne 0 ]; then
        echo "minimap2 failed" >&2
        exit 1
    fi

	return 0
}

run_samtools() {
	echo '-----------------------------------------------'
	echo 'Running samtools'
	echo '-----------------------------------------------'

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sam" -o "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sorted.bam"

	rm -f "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sorted.bam"

	if [ $? -ne 0 ]; then
        echo "samtools failed" >&2
        exit 1
    fi

	return 0
}

run_blobtools() {
	echo '-----------------------------------------------'
	echo 'Running blobtools'
	echo '-----------------------------------------------'

    conda activate "${BLOBTOOLS_CONDA}"

	if [ $? -ne 0 ]; then
        echo "Failed to activate blobtools conda environment" >&2
        exit 1
    fi

    ${BLOBTOOLS} create -i "${RUN_DIR}/assembly.fasta" -b "${RUN_DIR}/minimap2/${SAMPLE}_Raven.sorted.bam" \
		-t "${RUN_DIR}/contigs_vs_nt.blast" -o "${RUN_DIR}/blobtools" --db "${BLOBTOOLS_DB}"

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

	## Run Raven to assemble the reads
    run_raven

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
