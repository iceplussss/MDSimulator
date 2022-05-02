#!/bin/bash -l
#SBATCH --job-name=md
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

module load anaconda3/2021.5  
conda activate ase
python simulation.py 