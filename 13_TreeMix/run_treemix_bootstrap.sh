#!/usr/bin/env bash
# =============================================================================
# Run TreeMix bootstrap replicates for one fixed migration value, m=8.
# Each SLURM array task produces one independent bootstrap replicate.
# =============================================================================

#SBATCH --cluster=mindwell
#SBATCH --job-name=treemix_boot_m8
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/04_treemix_boot_m8_%A_rep%a.log
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-500%80

set -euo pipefail

# ── Environment ────────────────────────────────────────────────────────────────
module load cluster/wice/batch

source "/data/leuven/XXX/vscXXXXX/miniconda3/etc/profile.d/conda.sh"
conda activate treemix

# ── Working directory ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/treemix


# ── Directories and files ─────────────────────────────────────────────────────
INDIR="treemix_input"
OUTDIR="treemix_results"

INPUT="${INDIR}/treemix.frq.gz"
RUNDIR="${OUTDIR}/treemix_bootstrap_m8"

# ── TreeMix parameters ────────────────────────────────────────────────────────
ROOT_POP="CN_W1"

M=8
NBOOT=500
TREEMIX_K=50
GLOBAL_REARR=true

SEED_BASE=42000

# ── Bootstrap replicate handled by this array task ────────────────────────────
rep="${SLURM_ARRAY_TASK_ID}"
seed=$((SEED_BASE + rep))

mkdir -p "$RUNDIR"

outpfx="${RUNDIR}/treemix.m${M}.boot${rep}"

# ── Report settings ───────────────────────────────────────────────────────────
echo "Migration edges: $M"
echo "Bootstrap replicate: $rep of $NBOOT"
echo "Seed: $seed"
echo "Input frequencies: $INPUT"
echo "Output prefix: $outpfx"
echo "Root population: $ROOT_POP"
echo "TreeMix block size: $TREEMIX_K"
echo "Global rearrangements: $GLOBAL_REARR"

# ── Construct TreeMix command ─────────────────────────────────────────────────
cmd=(
    treemix
    -i "$INPUT"
    -root "$ROOT_POP"
    -o "$outpfx"
    -m "$M"
    -k "$TREEMIX_K"
    -seed "$seed"
    -bootstrap
)

[[ "$GLOBAL_REARR" == "true" ]] && cmd+=(-global)


# ── Run TreeMix ───────────────────────────────────────────────────────────────

echo "Running TreeMix bootstrap: m=$M replicate=$rep seed=$seed"

"${cmd[@]}" > "${outpfx}.log" 2>&1

echo "Bootstrap replicate $rep complete"
