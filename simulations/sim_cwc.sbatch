#!/bin/bash -l
#SBATCH --job-name=cwct4
#SBATCH --chdir="/beegfs/users/gmastrantonio/tokyo/"
#SBATCH --partition=reserved-gm
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000 
#SBATCH --output="res/%x_%j.out"
#SBATCH --error="res/%x_%j.err.out"



module load R/4.4.1-gfbf-2023b 


Rscript --vanilla  codes/cwc_simulations_large_d.R "$@"


