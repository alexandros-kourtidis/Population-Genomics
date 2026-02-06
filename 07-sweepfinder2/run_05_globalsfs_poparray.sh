#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=sf2_glsfs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-20		#20 populations

set -euo pipefail

# -------------------------------------------------------------------------------
# Create a global/empirical SFS from the filtered input files, per population.
# -------------------------------------------------------------------------------
# Environment setup
# -------------------------------------------------------------------------------
module load cluster/genius/amd
module load HTSlib/1.17-GCC-12.2.0
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh
conda activate py27
echo "================="

# -------------------------------------------------------------------------------
# Variables
# -------------------------------------------------------------------------------

cd /lustre1/scratch/363/vsc36396/sf2

scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)
populations=( A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM \
		B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 \
		C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 \
		D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050)

ID=$((SLURM_ARRAY_TASK_ID - 1))
pop="${populations[$ID]}"

IN_DIR="input_files/${pop}"
FILTER_DIR="${IN_DIR}/filtered"

header_source="${FILTER_DIR}/out.${pop}_${scaffolds[0]}.sweepfinder.filtered.input"
GLOBAL_ALLSCAFF="${FILTER_DIR}/${pop}_allscaffolds.sweepfinder.filtered.input"
GLOBAL_SFS="${FILTER_DIR}/${pop}_pol_globalsfs.input"

# -------------------------------------------------------------------------------
# Create the global/empirical SFS file once per population (if missing or empty)
# -------------------------------------------------------------------------------
#  Build combined input with header from unfiltered first scaffold and all filtered bodies
#  The tail -n +2 removes the header line from filtered inputs
# -------------------------------------------------------------------------------

echo "[${pop}] Creating global/empirical SFS files…"

{ head -n 1 "$header_source"
	tail -n +2 -q "${FILTER_DIR}"/out.${pop}_scaffold_*.sweepfinder.filtered.input
	} > "$GLOBAL_ALLSCAFF"

	# Run SweepFinder2 to make the global SFS
/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -f \
	"$GLOBAL_ALLSCAFF" \
	"$GLOBAL_SFS"

echo "[${pop}] Global/Empirical SFS file created: $GLOBAL_SFS"
