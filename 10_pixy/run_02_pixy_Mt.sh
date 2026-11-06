#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=pixy
#SBATCH --output=logs/pixy.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1             #1 mitochondrion

# -------------------------------------------------------------------------------
#  Running pixy for pi, theta w, and Tajima�s D for each population and
#  dxy, FST for each populaion pair.
#  Sources:
#	https://pixy.readthedocs.io/en/latest/index.html
#	https://doi.org/10.1111/1755-0998.13326
# -------------------------------------------------------------------------------
#  Environment variables
# -------------------------------------------------------------------------------
cd /lustre1/scratch/XXX/vscXXXXX/pixy
conda activate pixy

	# 1. Scaffold list
scaffolds=("Mt")
SCAFFOLD=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}

	# 2. Folders
IN_DIR="/lustre1/scratch/XXX/vscXXXXX/VCFs"
OUT_DIR="pixy_out_Mt"
mkdir -p ${OUT_DIR}

	# 3. File names & variables 
VCF="$IN_DIR/batch1234_Mt.invsnps.noMitoMiss_DP_MISS_haplo.vcf.gz"        
POPFILE="/user/leuven/XXX/vscXXXXX/DATA/scripts/10_pixy/infofile.txt"

WINDOW=15000        # bigger than the mt size = a single metric for the entire mitochondrion                             
STATS="pi fst dxy watterson_theta tajima_d"
fst_type="hudson"

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
     --output_prefix pixy_output_${SCAFFOLD}_${WINDOW} \
     --fst_type ${fst_type} \
     --fst_components \
     --tajima_components 

echo "pixy has run"
