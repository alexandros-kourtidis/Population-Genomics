#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=picmin_agri
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 	#10 scaffolds

# ========================================================================================
# PicMin SLURM Array
# ----------------------------------------------------------------------------------------
# This SLURM job array launches the PicMin analysis separately for each scaffold listed in the "scaffolds"
# Bash array below. The array index determines which scaffold is processed in each task.
# All parameters required by the R script are defined in this file and passed positionally to the R script
# using command-line arguments.
# ----------------------------------------------------------------------------------------
# IMPORTANT:
#  - Population names (POPS): 	a single, comma-separated string (no whitespace outside commas).
#  - Missing-levels (MISS): 	a single comma-separated string of integers.
#  - Input file path (INFILE): 	the global data file containing *all* scaffolds, since the R script
#     				will subset data internally based on the scaffold name.
#  - OUTPREFIX:			the first part of the output filename. The final output file will be:
#     				<OUT_DIR>/<OUTPREFIX><SCAFFOLD>.tsv
# ----------------------------------------------------------------------------------------
# Argument order passed to R is strictly important!!
# SCAFFOLD POPS MISS IN_DIR OUT_DIR INFILE OUTPREFIX WINDOW NUMREPS NSIMS
# ----------------------------------------------------------------------------------------
# All required parsing, validation, and output construction is performed inside the R script.
# ========================================================================================

set -euo pipefail

# 0. Environment
conda activate picmin

# 1. Scaffold list
scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
SCAFFOLD=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}

# 2. Populations
POPS="A1_AnOudin,A2_BuSN,A3_Mech,A4_PL15_YEL,A5_GenM,B1_BKN1,B2_OM2,B3_ZW,B4_OHZ,B5_DA2,C1_BlfN,C2_MO,C3_Ter1,C4_BW_48630,C5_BW_36962,D1_CBOO6,D2_LRV,D3_BKLE5,D4_BW_62256,D5_BW_22050"

# 3. Missing levels
MISS="10,11,12,13,14,15,16,17,18,19,20"

# 4. Dirs
IN_DIR="/lustre1/scratch/363/vsc36396/picmin/input"
OUT_DIR="/lustre1/scratch/363/vsc36396/picmin/output"

# 5. Input file
INFILE="${IN_DIR}/PicMin_input_meanLR_10kb_windows_woContam.txt"
OUTPREFIX="pops_all_"

# 6. Parameters
WINDOW=10000
NUMREPS=100000
NSIMS=40000

# ---- Run the R script ----

Rscript run_08_picmin_onescaf.R \
    "$SCAFFOLD" \
    "$POPS" \
    "$MISS" \
    "$IN_DIR" \
    "$OUT_DIR" \
    "$INFILE" \
    "$OUTPREFIX" \
    "$WINDOW" \
    "$NUMREPS" \
    "$NSIMS"
