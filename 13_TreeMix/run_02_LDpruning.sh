#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=ldpruning
#SBATCH --output=logs/ldpruning10kb.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# =================================================================================
#  LD pruning with PLINK 1.9 and turning into VCF for TreeMix downstream
# =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load PLINK/1.9
module load BCFtools/1.18-GCC-12.3.0

cd /lustre1/scratch/XXX/vscXXXXX/treemix

# ── Directories & Files ─────────────────────────────────────────────────────────
VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/treemix/VCFs"
OUT_DIR="plink_out"
mkdir -p "${OUT_DIR}"

VCF="${VCF_DIR}/batch1234CN.snps.MISS.vcf.gz"
PLINK_BASE="${OUT_DIR}/batch1234CN.snps.MISS.data"
OUT="${OUT_DIR}/batch1234CN.snps.MISS_plink10kb"

sum="${OUT_DIR}/batch1234CN.snps.MISS_plink10kb.summary.txt"

# ── Settings ───────────────────────────────────────────────────────
WINDOW_SIZE="10 kb"		# derived by the LD decay curved
STEP_SIZE=1				# should be 1 when window size is in kb
R2=0.1					# considered a strict filter, derived by the LD decay curve for all pops

# -------------------------------------------------------------------------------
#  Run the LD pruning 
# -------------------------------------------------------------------------------
echo "Running PLINK..."

echo "Converting VCF to PLINK binary..."
plink --vcf "${VCF}" \
	--double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --keep-allele-order \
	--make-bed \
	--out "${PLINK_BASE}"

echo "Calculating the correlation values..."
plink --bfile "${PLINK_BASE}" \
	--allow-extra-chr --allow-no-sex \
	--indep-pairwise \
		"${WINDOW_SIZE}" \
		"${STEP_SIZE}" \
		"${R2}" \
	--out "${OUT}"

echo "Creating the bed files..."
plink --bfile "${PLINK_BASE}" \
	--allow-extra-chr --allow-no-sex \
	--extract "${OUT}.prune.in" \
	--make-bed \
	--out "${OUT}_LDpruned"

echo "Converting PLINK binary back to vcf..."
plink --bfile "${OUT}_LDpruned" \
	 --allow-extra-chr --allow-no-sex --keep-allele-order \
	 --recode vcf-iid bgz \
	 --out "${OUT}"

echo "Indexing the vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "${OUT}.vcf.gz"

# Creating stats of the filtered vcf file
echo "creating stats txt file for ${OUT}.vcf.gz..."
bcftools stats -s - "${OUT}.vcf.gz" > "${sum}"

echo "Plink has run"