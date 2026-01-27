#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=sweepfinder2_A1_AnOudin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 	# 10 scaffolds

# -------------------------------------------------------------------------------
#	Runs SweepFinder2 on one population, on an array for all chromosomes
#   Note:
#	You need to run the first chromosome once, in order to create the global sfs
#	file, and then you can run all the rest chromosomes as an array.
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

IN_DIR="input_files"
OUT_DIR="output_files"

GLOBAL_SFS="$IN_DIR/${pop}_globalsfs.input"
GLOBAL_ALLSCAFF="$IN_DIR/${pop}_allscaffolds.sweepfinder.input"

mkdir -p "$OUT_DIR"

# -------------------------------------------------------------------------------
#	Create the global SFS file only if it does not already exist
# -------------------------------------------------------------------------------

if [[ ! -f "$GLOBAL_SFS" || ! -f "$GLOBAL_ALLSCAFF" ]]; then
    echo "Global SFS files not found — creating them now..."

{ head -n 1 $IN_DIR/out.${pop}_${scaffolds[0]}.sweepfinder.input; \
	tail -n +2 -q $IN_DIR/*.input; } > $IN_DIR/${pop}_allscaffolds.sweepfinder.input  
 
/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -f \
	"$GLOBAL_ALLSCAFF" \
	"$GLOBAL_SFS"

echo "Global SFS files created."
else
echo "Global SFS already exists — skipping creation."
fi

# -------------------------------------------------------------------------------
#	Run sweepfinder2 with global sfs
# -------------------------------------------------------------------------------

/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -lg 500 \
	$IN_DIR/out.${pop}_${scaf}.sweepfinder.input \
	$IN_DIR/${pop}_globalsfs.input \
	$OUT_DIR/${pop}_${scaf}_wglsfs.sf2.out 

echo "Sweepfinder2 has run for ${scaf} of ${pop}"
