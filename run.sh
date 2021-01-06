#!/bin/bash

#SBATCH --job-name=MAGITICS_kmers
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --ntasks-per-node=5
#SBATCH --array=10

module load miniconda
module load dsk
conda activate /home/ylucas/.conda/env/python3_yvan/
pip3 install xlrd
python3 -u magitics/main.py --len_kmers $SLURM_ARRAY_TASK_ID > run.txt 

