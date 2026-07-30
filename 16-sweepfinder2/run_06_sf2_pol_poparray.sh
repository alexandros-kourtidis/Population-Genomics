#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=sf2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-20 		#20 populations

set -euo pipefail

# -------------------------------------------------------------------------------
# Runs polarized SweepFinder2 per population, iterating over all scaffolds.
# Note:
# - sf2 doesn't support multithreading, so we loop sequentially over scaffolds.
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Environment setup
# -------------------------------------------------------------------------------

module load cluster/genius/amd
module load HTSlib/1.17-GCC-12.2.0
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh
conda activate py27
echo "================="

cd /lustre1/scratch/363/vsc36396/sf2

# -------------------------------------------------------------------------------
# Variables
# -------------------------------------------------------------------------------

	# Chromosome list (order matters for header source)
scaffolds=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10)

	# Population IDs 
populations=( \
	A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM \
	B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 \
	C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 \
	D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050\
)

# Resolve this array task's population (1-based indexing from SLURM_ARRAY_TASK_ID)
pop_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( pop_index < 0 || pop_index >= ${#populations[@]} )); then
    echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} out of range (1..${#populations[@]})"
    exit 1
fi
pop="${populations[$pop_index]}"

IN_DIR="input_files/${pop}"
FILTER_DIR="${IN_DIR}/filtered"
OUT_DIR="output_files/${pop}"

mkdir -p "$OUT_DIR"

GLOBAL_SFS="${FILTER_DIR}/${pop}_pol_globalsfs.input"

# -------------------------------------------------------------------------------
# Run SweepFinder2 with the global SFS in a 'for' loop for each scaffold
# -------------------------------------------------------------------------------

echo "Running SweepFinder2 for ${scaffolds[@]} in population ${pop}"

for scaf in "${scaffolds[@]}"; do
    FILTERED_FILE="${FILTER_DIR}/out.${pop}_${scaf}.sweepfinder.filtered.input"
    OUT_FILE="${OUT_DIR}/${pop}_${scaf}_wglsfs.sf2.out"

    if [[ ! -s "$FILTERED_FILE" ]]; then
        echo "[${pop} | ${scaf}] WARNING: Missing filtered input: $FILTERED_FILE — skipping."
        continue
    fi

    if [[ -s "$OUT_FILE" ]]; then
        echo "[${pop} | ${scaf}] Output exists, replacing: $OUT_FILE"
    fi

    echo "[${pop} | ${scaf}] Running SweepFinder2…"
    /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/07_sweepfinder2/SF2/SweepFinder2 -lg 500 \
        "$FILTERED_FILE" \
        "$GLOBAL_SFS" \
        "$OUT_FILE"

    echo "[${pop} | ${scaf}] Done."
done

echo "[${pop}] Completed all scaffolds."
