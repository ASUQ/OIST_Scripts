#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2
#SBATCH --mem=100G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --job-name=busco_multigene_phylogeny
#SBATCH --output=/path/to/dir/busco_multigene_phylogeny.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user.name@oist.jp

#-------------Softwares--------------#
# make sure you installed miniconda3
source '~/.bashrc'

# conda create -n busco_phylogeny \
#   Biopython, tqdm, mafft, trimAl, AMAS, IQ-TREE

CONDA_ENV='busco_phylogeny'

SCRIPT='./Busco_multigene_tree.py'
#------------------------------------#


INPUT_DIR="<INPUT_DIR>"
OUTPUT_DIR="<OUTPUT_DIR>"

# Example
conda activate "${CONDA_ENV}"
"${SCRIPT}" all --input_dir "${INPUT_DIR}" --out_dir "${OUTPUT_DIR}"  \
    --cores 16 --fraction '0.80,0.85,0.90,0.95,1.00' --verbose \
    --mafft "--localpair --maxiterate 1000 --thread 16"

conda deactivate
