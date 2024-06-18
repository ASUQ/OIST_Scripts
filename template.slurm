#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --job-name=MMDD_Job
#SBATCH --output=/path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-name@oist.jp

### Slurm information: https://groups.oist.jp/scs/use-slurm
### Save as /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm

#===============================================================================================#

# Description:
# This script is a template for running job with Deigo (computer cluster at OIST).
#
# Workflow:
#
# Expected Inputs:
#
# Required Softwares:
#
# Usage:

#===============================================================================================#

set -euxo pipefail

#===============================================================================================#

## Edit this part before running the script

#-------------Softwares--------------#
# Conda env paths
MINICONDA_DIR='Path to miniconda directory'

# Software paths

# Database paths

#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='<Path to directory for running job>'
#------------------------------------#

#-------------Parameters-------------#

#------------------------------------#


#===============================================================================================#

# DON'T MODIFY THE SCRIPT BELOW!!

activate_conda() {
    source "${MINICONDA_DIR}/etc/profile.d/conda.sh"
    if ! command -v conda &> /dev/null; then
        echo "conda could not be found"
        exit 1
    fi
}

main() {
    #＃ Add conda to your PATH
	activate_conda

	## Create running directory if doesn't exist
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}" || exit

	# module installation

	# Run functions

    echo 'JOB COMPLETED'
}


# Execute main function
main

### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
