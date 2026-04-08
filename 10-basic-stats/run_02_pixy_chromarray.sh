#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=pixy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-10 	#10 scaffolds

# -------------------------------------------------------------------------------
#  Running pixy for pi, theta w, and Tajima’s D for each population and
#  dxy, FST for each populaion pair.
#  Sources:
#	https://pixy.readthedocs.io/en/latest/index.html
#	https://doi.org/10.1111/1755-0998.13326
# -------------------------------------------------------------------------------
#  Environment variables
# -------------------------------------------------------------------------------
cd /lustre1/scratch/363/vsc36396/fst
conda activate pixy

	# 1. Scaffold list
scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
SCAFFOLD=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}

	# 2. Folders
IN_DIR="filtered_vcfs"
OUT_DIR="pixy_out"
mkdir -p ${OUT_DIR}

	# 3. File names & variables 
VCF="$IN_DIR/batch1234.filtered_w_inv.vcf.gz"        
POPFILE="population_map.txt"

WINDOW=10000                              
STATS="pi fst dxy watterson_theta tajima_d"

# -------------------------------------------------------------------------------
# 1. Run pixy
# 	All possible arguments for stats: pi fst dxy watterson_theta tajima_d
# -------------------------------------------------------------------------------
echo "Running pixy for ${STATS} on ${SCAFFOLD} every ${WINDOW}bp ..."

pixy --stats ${STATS} \
     --vcf ${VCF} \
     --populations ${POPFILE} \
     --chromosomes ${SCAFFOLD} \
     --window_size ${WINDOW} \
     --n_cores ${SLURM_CPUS_PER_TASK} \
     --output_folder ${OUT_DIR} \
     --output_prefix pixy_output_${SCAFFOLD}_${WINDOW}

echo "pixy has run"
