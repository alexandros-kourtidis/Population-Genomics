#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=ldpruning
#SBATCH --output=logs/ldpruning20kb.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# =================================================================================
#  LD pruning with PLINK 1.9
# =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load PLINK/1.9
  
cd /lustre1/scratch/XXX/vscXXXXX/ld_pruning

# ── Directories & Files ─────────────────────────────────────────────────────────
VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/VCFs"
OUT_DIR="plink_out"
mkdir -p "${OUT_DIR}"

VCF="${VCF_DIR}/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz"
OUT="${OUT_DIR}/batch1234_plink20kb"

# ── Settings ───────────────────────────────────────────────────────
WINDOW_SIZE="20 kb"		# derived by the LD decay curved
STEP_SIZE=1				# should be 1 when window size is in kb
R2=0.1					# considered a strict filter, derived by the LD decay curve for all pops

# -------------------------------------------------------------------------------
#  Run the LD pruning 
# -------------------------------------------------------------------------------
echo "Running PLINK..."

echo "Converting VCF to PLINK binary..."
plink --vcf "${VCF}" \
	--double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
	--make-bed \
	--out data

echo "Calculating the correlation values..."
plink --bfile data \
	--allow-extra-chr --allow-no-sex \
	--indep-pairwise \
		"${WINDOW_SIZE}" \
		"${STEP_SIZE}" \
		"${R2}" \
	--out "${OUT}"

echo "Creating the bed files..."
plink --bfile data \
	--allow-extra-chr --allow-no-sex \
	--extract "${OUT}.prune.in" \
	--make-bed \
	--out "${OUT}_LDpruned"

echo "Plink has run"