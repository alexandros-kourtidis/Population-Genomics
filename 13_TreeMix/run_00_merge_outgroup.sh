#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=vcf_merge
#SBATCH --output=logs/vcf_merge.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
 
# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/treemix

OUT_DIR="VCFs"
mkdir -p "${OUT_DIR}"

VCF1="/lustre1/scratch/XXX/vscXXXXX/VCFs/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz"
VCF2="/lustre1/scratch/XXX/vscXXXXX/CN_W1/CN_DmagnaLRV01_merged1_10.vcf.gz"
OUTVCF="${OUT_DIR}/batch1234CN_DmagnaLRV01_merged1_10.vcf.gz"

sum="${OUT_DIR}/batch1234CN_DmagnaLRV01_merged1_10.summary.txt"

# -------------------------------------------------------------------------------
#  Merging 
# -------------------------------------------------------------------------------

# Merge the VCFs of different samples ( into a single VCF file.
echo "merging ${VCF1} and ${VCF2} ..."
bcftools merge \
	--threads ${SLURM_CPUS_PER_TASK} \
	"${VCF1}" \
	"${VCF2}" \
	-o "${OUTVCF}" \
	-O z

echo "indexing the merged vcf..."
bcftools index -t "$OUTVCF"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $OUTVCF..."
bcftools stats -s - "$OUTVCF" > "$sum"

echo "merging is done"