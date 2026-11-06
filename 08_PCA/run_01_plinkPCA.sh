#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=pca
#SBATCH --output=logs/pca.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 

# =================================================================================
#  PCA with PLINK 1.9
# =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load PLINK/1.9

cd /lustre1/scratch/XXX/vscXXXXX/pca

# ── Directories & Files ─────────────────────────────────────────────────────────
IN_DIR="/lustre1/scratch/XXX/vscXXXXX/ld_pruning/plink_out"
OUT_DIR="plink_out"
mkdir -p "${OUT_DIR}"

IN_BED="${IN_DIR}/batch1234_plink20kb_LDpruned"
OUT_PCA="${OUT_DIR}/pruned_20kb0.1"

# ── Settings ───────────────────────────────────────────────────────
PCs=50

# -------------------------------------------------------------------------------
#  Run the pca
# -------------------------------------------------------------------------------
echo "Running the PCA in PLINK..."

plink --bfile "${IN_BED}" \
    --allow-extra-chr --allow-no-sex \
    --pca ${PCs}  \
    --out "${OUT_PCA}" 

echo "plink pca is done"