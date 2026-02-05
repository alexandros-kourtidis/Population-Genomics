#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=sf2_A1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-10 	# !!! ALWAYS run array=1 first and then the rest - to ensure global SFS creation !!! #

# -------------------------------------------------------------------------------
#	Runs polarised SweepFinder2 on one population, on an array for all chromosomes
#   Notes:
#	- You need to run the first chromosome once, in order to create the global sfs
#	file, and then you can run all the rest chromosomes as an array.
#	- sf2 doesn't support multithreading.
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
#	Environment setup
# -------------------------------------------------------------------------------

module load cluster/genius/amd
module load HTSlib/1.17-GCC-12.2.0
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh 
conda activate py27
echo "================="

cd /lustre1/scratch/363/vsc36396/sf2

# -------------------------------------------------------------------------------
#	Variables assignment
# -------------------------------------------------------------------------------

scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
scaf=${scaffolds[$((SLURM_ARRAY_TASK_ID - 1))]}
pop="A1_AnOudin"

IN_DIR="input_files/${pop}"
FILTER_DIR="${IN_DIR}/filtered"
OUT_DIR="output_files/${pop}"

mkdir -p "$OUT_DIR"

FILTERED_FILE="$FILTER_DIR/out.${pop}_${scaf}.sweepfinder.filtered.input"
OUT_FILE="$OUT_DIR/${pop}_${scaf}_wglsfs.sf2.out"

GLOBAL_ALLSCAFF="$IN_DIR/${pop}_allscaffolds.sweepfinder.input"
GLOBAL_SFS="$IN_DIR/${pop}_pol_globalsfs.input"

# -------------------------------------------------------------------------------
#	Create the global/empirical SFS file
# -------------------------------------------------------------------------------

echo "Creating global/empirical SFS file"

{ head -n 1 $IN_DIR/out.${pop}_${scaffolds[0]}.sweepfinder.input; \
	tail -n +2 -q "$FILTER_DIR"/out.${pop}_*.sweepfinder.filtered.input; } > "$GLOBAL_ALLSCAFF"  
 
/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -f \
	"$GLOBAL_ALLSCAFF" \
	"$GLOBAL_SFS"

echo "Global/Empirical SFS file created."

# -------------------------------------------------------------------------------
#	Run sweepfinder2 with global/empirical sfs
# -------------------------------------------------------------------------------

/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -lg 500 \
	"$FILTERED_FILE" \
	"$GLOBAL_SFS" \
	"$OUT_FILE"

echo "Sweepfinder2 has run for ${scaf} of ${pop}"
