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
# 1. Find read overlaps with minimap2
# 2. Assemble long read with miniasm
# 3. Annotate scaffolds taxonomically with megablast
# 4. Map reads onto assembly with minimap2 and samtools
# 5. Summarize results with Blobtools
#
# Expected Inputs:
# - Pacbio or Oxford Nanopore long read    # You need to adjust the Flye option based on the read type
#
# Required Softwares:
# - miniasm
# - any2fasta
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

## Edit this part before running the script

#-------------Softwares--------------#
# Conda env paths
MINICONDA_DIR='Path to miniconda directory'
BLOBTOOLS_CONDA='Path to blobtools conda environment'
MINIASM_CONDA='Path to miniasm conda environment'	# I expect any2fasta is installed in same conda env.

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
}

activate_conda() {
    source "${MINICONDA_DIR}/etc/profile.d/conda.sh"
    if ! command -v conda &> /dev/null; then
        echo "conda could not be found"
        exit 1
    fi
}

run_minimap2_overlap() {
	echo '-----------------------------------------------'
	echo 'Running minimap2 to overlap reads'
	echo '-----------------------------------------------'

	# Adjust the read option based on your data
	# - ava-pb/ava-ont - PacBio/Nanopore read overlap

	${MINIMAP2} -x ava-ont -t $((THREADS)) "${long_read}" "${long_read}" | gzip > "${RUN_DIR}/overlaps.paf.gz"

    if [ $? -ne 0 ]; then
        echo "minimap2 overlap failed" >&2
        exit 1
    fi
}

run_miniasm() {
    echo '-----------------------------------------------'
	echo 'Running miniasm'
	echo '-----------------------------------------------'

	conda activate "${MINIASM_CONDA}"

	if [ $? -ne 0 ]; then
        echo "Failed to activate miniasm conda environment" >&2
        exit 1
    fi

	miniasm -f "${long_read}" "${RUN_DIR}/overlaps.paf.gz" > "${RUN_DIR}/assembly.gfa"

    if [ $? -ne 0 ]; then
        echo "miniasm failed" >&2
        exit 1
    fi

	any2fasta "${RUN_DIR}/assembly.gfa" > "${RUN_DIR}/assembly.fasta"

    if [ $? -ne 0 ]; then
        echo "any2fasta failed" >&2
        exit 1
    fi

	conda deactivate
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
}

run_minimap2() {
	echo '-----------------------------------------------'
	echo 'Running minimap2'
	echo '-----------------------------------------------'

	mkdir -p minimap2

	# Adjust the read option based on your data
	# - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
	# - map-hifi - PacBio HiFi reads vs reference mapping

	${MINIMAP2} -a -x map-ont -t $((THREADS)) "${RUN_DIR}/assembly.fasta" "${long_read}" > "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sam"

    if [ $? -ne 0 ]; then
        echo "minimap2 failed" >&2
        exit 1
    fi
}

run_samtools() {
	echo '-----------------------------------------------'
	echo 'Running samtools'
	echo '-----------------------------------------------'

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sam" -o "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sorted.bam"

	rm -f "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sorted.bam"

	if [ $? -ne 0 ]; then
        echo "samtools failed" >&2
        exit 1
    fi
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

    ${BLOBTOOLS} create -i "${RUN_DIR}/assembly.fasta" -b "${RUN_DIR}/minimap2/${SAMPLE}_miniasm.sorted.bam" \
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
}

main() {
    #＃ Add conda to your PATH
	activate_conda

	## Create running directory if doesn't exist
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}" || exit

	## Set software paths
	set_software_paths

	## Run minimap2 to find read overlaps
	run_minimap2_overlap

	## Run miniasm to assemble the reads
    run_miniasm

	## Run megablast to taxonomically annotate the scaffolds
    run_blastn

	## Run minimap2 to map the read to the reference
	run_minimap2

	## Run samtools to convert sam file to sorted bam file and index
	run_samtools

	## Run blobtools to summarize the blastn result
    run_blobtools

    echo 'JOB COMPLETED'
}


# Execute main function
main

### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
