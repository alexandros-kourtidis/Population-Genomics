#!/usr/bin/env bash
# =============================================================================
# LD pruning & PCA with PLINK
# (PLINK2 uses a different file format and is stuck in 'beta' version)
# =============================================================================
#SBATCH --cluster=mindwell
#SBATCH --job-name=filter_treemix
#SBATCH --output=logs/filter_MISS.%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH -A lp_svbelleghem

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/treemix

module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/treemix/VCFs"

IN_VCF="${VCF_DIR}/batch1234CN_DmagnaLRV01_merged1_10.vcf.gz"
FILT_VCF="${VCF_DIR}/batch1234CN.snps.MISS.vcf.gz"

sum="${VCF_DIR}/batch1234CN.snps.MISS.summary.txt"

# ── Settings ───────────────────────────────────────────────────────
MAX_MISSING=0.10        # maximum fraction of missing genotypes per site

# -------------------------------------------------------------------------------
# Filter with bcftools 
#    The vcf containing our dataset snps is already filtered (see 04_filtering: batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz)
#    I have merged the filtered snps with the unfiltered CN_W1 and then I filter for missingness to remove all sites not present in the filtered snp file. 
# -------------------------------------------------------------------------------

echo "filtering for missingness ..."
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 \
    -M2 \
    -v snps \
    -i "F_MISSING<=$MAX_MISSING" \
    -Oz \
    -o "$FILT_VCF" \
    "$IN_VCF" 

# Index the filtered VCF
echo "indexing the filtered vcf..."
bcftools index -t "$FILT_VCF"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $FILT_VCF..."
bcftools stats -s - "$FILT_VCF" > "$sum"

echo "filtering is done"
