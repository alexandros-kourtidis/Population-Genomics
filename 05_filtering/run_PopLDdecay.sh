#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=popLDdecay
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 	#20 populations

# -------------------------------------------------------------------------------
#  Running PopLDdecay
# -------------------------------------------------------------------------------
#  Environment variables
# -------------------------------------------------------------------------------
cd /lustre1/scratch/363/vsc36396/ld_decay
conda activate PopLDdecay

	# 1. Populations list
pops=(A1_AnOudin A2_BuSN A3_Mech A4_Pl15YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW48630 C5_BW36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW62256 D5_BW22050)
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}

	# 2. Folders
IN_DIR=".."
POP_DIR="pops_lists"
OUT_DIR="ld_decay_out"
mkdir -p ${OUT_DIR}

	# 3. File names & variables 
VCF="$IN_DIR/batch1234_DmagnaLRV01_merged1_10.filtered.vcf.gz"        
POPFILE="${POP_DIR}/${pop}_samples.txt"

# -------------------------------------------------------------------------------
# 1. Run PopLDdecay 
# -------------------------------------------------------------------------------

echo "Running PopLDdecay for population ${pop}..."

PopLDdecay \
    -InVCF ${VCF} \
    -OutStat ${OUT_DIR}/${pop}.LDdecay \
    -SubPop ${POPFILE} \
    -MaxDist 500000 \
    -MAF 0.05 \
    -Het 0.88 \
    -Miss 0.25

# -------------------------------------------------------------------------------
# 2. Plot LD decay 
# -------------------------------------------------------------------------------

echo "Plotting LD decay for population ${pop}..."

Plot_OnePop.pl \
	-inFile ${OUT_DIR}/${pop}.LDdecay.stat.gz \
	-binSize 5000 \
	-output ${OUT_DIR}/${pop}.LDdecay.png

echo "PopLDdecay has run"
