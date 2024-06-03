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

#-------------Softwares--------------#

#------------------------------------#


# Define file and directory paths
RUN_DIR='RUN_DIR'

# Create running directory if doesn't exist
mkdir -p $RUN_DIR
cd $RUN_DIR || exit

### Add your code here


### Run with $ sbatch /path/to/personal/directory/Jobs/YYYYMM/MMDD_Job.slurm
### To check the running jobs, run $ squeue
