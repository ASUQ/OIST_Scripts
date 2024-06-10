#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --job-name=MMDD_Job
#SBATCH --output=/path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.name@oist.jp

### Slurm information: https://groups.oist.jp/scs/use-slurm
### Save as /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm

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

main() {
    #＃ Add conda to your PATH
	source "${MINICONDA_DIR}/etc/profile.d/conda.sh"

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
