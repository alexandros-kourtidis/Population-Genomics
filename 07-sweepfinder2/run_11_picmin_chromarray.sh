#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=picmin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=10 	#10 scaffolds

set -euo pipefail

conda activate picmin

	# Define scaffold
scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
scaff=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}

	# Run the R script for picmin by scaffold
Rscript run_08_picmin_1scaf.R $scaff
