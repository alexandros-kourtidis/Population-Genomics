#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=thining
#SBATCH --output=logs/thining10kb.%j.log
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
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/hierfstat

VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/VCFs"
OUT_DIR="plink_out"
mkdir -p "${OUT_DIR}"

IN_VCF="${VCF_DIR}/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz"
OUT="${OUT_DIR}/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss_thin10kb"

# ── Settings ───────────────────────────────────────────────────────
WINDOW_SIZE=10000		# derived by the LD decay curved

# -------------------------------------------------------------------------------
#  Run the thining 
# -------------------------------------------------------------------------------
echo "Running PLINK..."

echo "Converting VCF to PLINK binary and thining by ${WINDOW_SIZE} bp..."
plink --vcf "${IN_VCF}" \
	 --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
	 --bp-space ${WINDOW_SIZE} \
	 --make-bed \
	 --out "${OUT}"

echo "Converting PLINK binary back to vcf..."
plink --bfile "${OUT}" \
	 --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
	 --recode vcf-iid bgz\
	 --out "${OUT}"

echo "Indexing the vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "${OUT}.vcf.gz"

echo "Plink thining has by ${WINDOW_SIZE} bp run"