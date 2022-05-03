#!/bin/bash -l
#SBATCH --job-name=md
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load anaconda3/2021.5  
conda activate ase
python benchmark.py 