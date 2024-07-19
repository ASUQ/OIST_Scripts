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
# This script is for reassembling bacteria genome by extracting reads from metagenome assembly with Deigo (computer cluster at OIST).
#
# Workflow:
# 1. Extract short reads of bacteria with bowtie2
# 2. Extract long reads of bacteria with minimap2
# 3. Assemble extracted reads with Unicycler
# 4. Annotate scaffolds taxonomically with megablast
# 5. Map reads onto assembly with bowtie2 and samtools
# 6. Summarize results with Blobtools
#
# Expected Inputs:
# - Illumina pair-end reads (quality controlled)
# - Pacbio or Oxford Nanopore long read
#
# Required Softwares:
# - Unicycler
# - SPADES
# - Standalone BLAST
# - bowtie2
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

# Software directory paths
BLAST='<BLAST+ directory>'
BLOBTOOLS='<Blobtools directory>'
BOWTIE2='<bowtie2 directory>'
MINIMAP2='<minimap2 directory>'
RACON='<Racon directory>'
SAMTOOLS='<samtools directory>'
SPADES='<SPAdes directory>'
UNICYCLER='<Unicycler directory>'

# Database paths
BLOBTOOLS_DB='Path to blobtools nodesDB'
BLAST_NT_DB='Path to nt database from ncbi'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'

SAMPLE='<sample name>'

reference='<reference>.fasta'

read_f='<Illumina forward read>.fastq.gz'
read_r='<Illumina reverse read>.fastq.gz'
long_read='<Long read>.fastq.gz'	# Adjust the minimap2 option based on your read types
#------------------------------------#

#-------------Parameters-------------#
THREADS='<number of CPUs>'	# More than 16 is recommended
MEMORY='<size of RAM>'  # More than 100 GB is recommended
#------------------------------------#


#===============================================================================================#

# DON'T MODIFY THE SCRIPT BELOW!!

set_software_paths() {
	# blast
	BLASTN="${BLAST}/bin/blastn"
	MAKEBLASTDB="${BLAST}/bin/makeblastdb"
	TBLASTN="${BLAST}/bin/tblastn"

	# blobtools
	BLOBTOOLS="${BLOBTOOLS}/blobtools"

	# bowtie
	BOWTIE2_INDEX="${BOWTIE2}/bowtie2-build"
	BOWTIE2="${BOWTIE2}/bowtie2"

	# minimap2
	MINIMAP2="${MINIMAP2}/minimap2"

	# racon
	RACON="${RACON}/build/bin/racon"

	# samtools
	SAMTOOLS="${SAMTOOLS}/samtools"

	# spades
	SPADES="${SPADES}/bin/spades.py"

	# unicycler
	UNICYCLER="${UNICYCLER}/unicycler-runner.py"

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

run_bowtie2_extract_reads() {
	echo '-----------------------------------------------'
	echo 'Running bowtie2 to extract reads'
	echo "$(date)"
	echo '-----------------------------------------------'

	mkdir -p extracted_reads/short_reads

	${BOWTIE2_INDEX} -f -q --threads $((THREADS)) "${reference}" "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads"
	${BOWTIE2} -q --reorder --threads $((THREADS)) --time --met-stderr --met 10 \
		-x "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads" -1 "${read_f}" -2 "${read_r}" \
		--al-conc-gz "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.fastq.gz" > /dev/null

    if [ $? -ne 0 ]; then
        echo "bowtie2 failed" >&2
        exit 1
    fi

	rm -f "${RUN_DIR}/bowtie2/${SAMPLE}_reads.*.bt2"

	mv "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.fastq.1.gz" "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.1.fastq.gz"
	mv "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.fastq.2.gz" "${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.2.fastq.gz"

	read_ex_f="${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.1.fastq.gz"
	read_ex_r="${RUN_DIR}/extracted_reads/short_reads/${SAMPLE}_reads.2.fastq.gz"

	return 0
}

run_minimap2_extract_reads() {
	echo '-----------------------------------------------'
	echo 'Running minimap2 to extract reads'
	echo "$(date)"
	echo '-----------------------------------------------'

	mkdir -p extracted_reads/long_reads

	# Adjust the read option based on your data
	# - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
	# - map-hifi - PacBio HiFi reads vs reference mapping

	${MINIMAP2} -a -x map-ont -t $((THREADS)) "${reference}" "${long_read}" > "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sam"

    if [ $? -ne 0 ]; then
        echo "minimap2 failed" >&2
        exit 1
    fi

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sam" -o "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sorted.bam"

	rm -f "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sorted.bam"

	# Extract mapped read from bam file
	${SAMTOOLS} fastq -@ $((THREADS-1)) -F 4 "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sorted.bam" | gzip > "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.fastq.gz"

    if [ $? -ne 0 ]; then
        echo "samtools2 failed" >&2
        exit 1
    fi

	rm -f "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sam"
	rm -f "${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.sorted.bam"

	long_ex_read="${RUN_DIR}/extracted_reads/long_reads/${SAMPLE}_reads.fastq.gz"

	return 0
}

run_unicycler() {
    echo '-----------------------------------------------'
	echo 'Running Unicycler'
	echo "$(date)"
	echo '-----------------------------------------------'

	# module installation
	module purge
	module load python/3.11.4	# Default python version on Deigo is v3.6.8, v3.8+ is needed for SPAdes v4.0+

	${UNICYCLER} --short1 "${read_ex_f}" --short2 "${read_ex_r}" --long "${long_ex_read}" --out "${RUN_DIR}" \
		--spades_path "${SPADES}" --makeblastdb_path "${MAKEBLASTDB}" --tblastn_path "${TBLASTN}" --racon_path "${RACON}" \
		--threads $((THREADS)) \--spades_options "--memory $((MEMORY))"

    if [ $? -ne 0 ]; then
        echo "Unicycler failed" >&2
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

    ${BLASTN} -task megablast -query "${RUN_DIR}/assembly.fasta" \
		-db "${BLAST_NT_DB}" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
		-culling_limit 5 -num_threads $((THREADS)) -evalue 1e-25 -max_target_seqs 5 > "${RUN_DIR}/assembly_vs_nt.blast"

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

	mkdir -p bowtie2

	${BOWTIE2_INDEX} -f -q --threads $((THREADS)) "${RUN_DIR}/assembly.fasta" "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler"
	${BOWTIE2} -q --reorder --threads $((THREADS)) --time --met-stderr --met 10 \
		-x "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler" -1 "${read_f}" -2 "${read_r}" > "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler".sam

    if [ $? -ne 0 ]; then
        echo "bowtie2 failed" >&2
        exit 1
    fi

	rm -f "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.*.bt2"

	return 0
}

run_samtools() {
	echo '-----------------------------------------------'
	echo 'Running samtools'
	echo "$(date)"
	echo '-----------------------------------------------'

	${SAMTOOLS} sort --output-fmt BAM --threads $((THREADS-1)) "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.sam" -o "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.sorted.bam"

	rm -f "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.sam"

	${SAMTOOLS} index -b -@ $((THREADS-1)) "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.sorted.bam"

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

    ${BLOBTOOLS} create -i "${RUN_DIR}/${SAMPLE}_Unicycler.contigs.fasta" -b "${RUN_DIR}/bowtie2/${SAMPLE}_Unicycler.sorted.bam" \
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

	## Run bowtie2 to extract reads
	run_bowtie2_extract_reads

	## Run minimap2 to extract reads
	run_minimap2_extract_reads

	## Run Unicycler to assemble the reads
    run_unicycler

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
