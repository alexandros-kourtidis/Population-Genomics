#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=admixture
#SBATCH --output=logs/admixture.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-30

# =================================================================================
#  ADMIXTURE
# ---------------------------------------------------------------------------------
# 	on SNPs filtered for quality, missingness, maf, and linkage prunned 
# ---------------------------------------------------------------------------------
# BEFORE RUNNING:
# 	ADMIXTURE does not accept chromosome names that are not human chromosomes. Will just exchange the first column by 0.
#	awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
#	mv $FILE.bim.tmp $FILE.bim
# =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
conda activate admixture
  
cd /lustre1/scratch/XXX/vscXXXXX/admixture

# ── Directories & Files ─────────────────────────────────────────────────────────
IN_DIR="/lustre1/scratch/XXX/vscXXXXX/ld_pruning/plink_out"
OUT_DIR="admixture_out"
mkdir -p "${OUT_DIR}"

IN="${IN_DIR}/batch1234_plink20kb_LDpruned"
OUT="${OUT_DIR}/pruned_20kb0.1"

# ── Settings ─────────────────────────────────────────────────────────
K=$SLURM_ARRAY_TASK_ID

# ── Run ADMIXTURE ─────────────────────────────────────────────────────────
echo "Running ADMIXTURE for K=${K} ..."
admixture --cv "${IN}.bed" \
	${K} \
	-j5 \
	> "${OUT}_log${K}.out"

echo "done"