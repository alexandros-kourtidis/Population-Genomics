#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=sweepfinder2_A1_AnOudin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-10 	# 10 scaffolds
 
# Load the programs we will use
module load cluster/genius/amd
module load HTSlib/1.17-GCC-12.2.0
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh 
conda activate py27

echo "================="

cd /lustre1/scratch/363/vsc36396/sweepfinder

	# The scaffolds of your genome and the population name
scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
scaf=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}
pop="A1_AnOudin"

	# Directories assignment
IN_DIR="input_files"
OUT_DIR="output_files"
mkdir -p "$OUT_DIR"

	# Create the global sfs file, only need to be done once!
# { head -n 1 $IN_DIR/${pop}_${scaffolds[1]}.sweepfinder.input; tail -n +2 -q $IN_DIR/*.input; } > $IN_DIR/${pop}_allscaffolds.sweepfinder.input  
# /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -f \
#	$IN_DIR/${pop}_allscaffolds.sweepfinder.input \
#	$IN_DIR/${pop}_globalsfs.input

	# Run sweepfinder2 with global sfs
/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -lg 500 \
	$IN_DIR/${pop}_${scaf}.sweepfinder.input \
	$IN_DIR/${pop}_globalsfs.input \
	$OUT_DIR/${pop}_${scaf}_wglsfs.sf2.out 

echo "Sweepfinder2 has run for ${scaf} of ${pop}"
