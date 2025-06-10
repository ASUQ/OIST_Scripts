#!/bin/bash
#SBATCH --partition short
#SBATCH --time 0-1
#SBATCH --mem=50G
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
# This script is to analyse taxonomical composition of metagenome using phyloflash with Deigo (computer cluster at OIST).
#
# Workflow:
# 1. Run Phyloflash
#
# Expected Inputs:
# - Illumina pair-end reads (quality controlled)
#
# Required Softwares:
# - Phyloflash
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
PHYLOFLASH_CONDA='Path to blobtools conda environment'

# Software paths

# Database paths
PHYLOFLASH_DB='Path to phyloflash DB'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'

read_f='<Illumina forward read>.fastq.gz'
read_r='<Illumina reverse read>.fastq.gz'

SAMPLE='sample_name'
#------------------------------------#

#-------------Parameters-------------#
THREADS="$SLURM_CPUS_PER_TASK"
# MEMORY="$SLURM_MEM_PER_NODE"
#------------------------------------#


#===============================================================================================#

# DON'T MODIFY THE SCRIPT BELOW!!

set_software_paths() {
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

run_phyloflash() {
    echo '-----------------------------------------------'
	echo 'Running Phyloflash'
	echo "$(date)"
	echo '-----------------------------------------------'

    conda activate "${PHYLOFLASH_CONDA}"

    phyloFlash.pl -lib "${SAMPLE}" -read1 "${read_f}" -read2 "${read_r}" \
        -dbhome "${PHYLOFLASH_DB}" -CPUs $((THREADS)) -everything > "${RUN_DIR}/phyloflash.output"

    conda deactivate

    if [ $? -ne 0 ]; then
        echo "Phyloflash failed" >&2
        exit 1
    fi


	return 0
}

main() {\
    ## Add conda to your PATH
	activate_conda

	## Create running directory if doesn't exist
    mkdir -p "${RUN_DIR}" && cd "${RUN_DIR}"

	## Set software paths
	set_software_paths

	## Run Phyloflash
    run_phyloflash

    echo 'JOB COMPLETED'

	return 0
}


# Execute main function
main

### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
