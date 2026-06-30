#!/usr/bin/env bash
# =============================================================================
# Run TreeMix for a single m value (0..MAX_M), submitted as a SLURM job array.
# Each array task handles one m, running REPLICATES sequentially.
# =============================================================================
#SBATCH --cluster=mindwell
#SBATCH --job-name=treemix_run
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/04_run_treemix_%A_m%a.log
#SBATCH -A lp_svbelleghem
#SBATCH --array=0    # amount of migration edges m

set -euo pipefail

echo "Loading cluster modules..."
module load cluster/wice/batch

source "/data/leuven/363/vsc36396/miniconda3/etc/profile.d/conda.sh"
conda activate treemix

# =============================================================================
# k=20 (TREEMIX_K=block grouping) because it is pre-pruned dataset
# replicates =1 for initial testing, then higher for (e.g.10) to test the topology of preferred m
# =============================================================================
cd /lustre1/scratch/363/vsc36396/treemix

INDIR="treemix_input"
OUTDIR="treemix_results"
INPUT="$INDIR/treemix.frq.gz"
RUNDIR="$OUTDIR/treemix_runs"

ROOT_POP="CN_W1"
REPLICATES=1
TREEMIX_K=20
GLOBAL_REARR=true

m="${SLURM_ARRAY_TASK_ID}"

mkdir -p "$RUNDIR/m${m}"

echo "================="
echo "m: $m"
echo "Input frequencies: $INPUT"
echo "Output dir: $RUNDIR/m${m}"
echo "Replicates: $REPLICATES"

SEED_BASE=42

for rep in $(seq 1 "$REPLICATES"); do
    seed=$((SEED_BASE + m * 1000 + rep))
    outpfx="$RUNDIR/m${m}/rep${rep}/treemix.m${m}.rep${rep}"
    mkdir -p "$(dirname "$outpfx")"

    cmd=(
        treemix
        -i "$INPUT"
        -root "$ROOT_POP"
        -o "$outpfx"
        -m "$m"
        -k "$TREEMIX_K"
        -seed "$seed"
    )

    [[ "$GLOBAL_REARR" == "true" ]] && cmd+=(-global)

    echo "Running m=$m rep=$rep seed=$seed"
    "${cmd[@]}" &> "${outpfx}.log"
done

echo "m=$m complete (all $REPLICATES replicates)"