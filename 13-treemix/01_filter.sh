#!/usr/bin/env bash
# =============================================================================
# LD pruning & PCA with PLINK
# (PLINK2 uses a different file format and is stuck in 'beta' version)
# =============================================================================
#SBATCH --cluster=mindwell
#SBATCH --job-name=filter_treemix
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/01_filter.%j.log
#SBATCH -A lp_svbelleghem

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────
cd /lustre1/scratch/363/vsc36396/treemix

module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
VCFDIR="/lustre1/scratch/363/vsc36396/treemix/vcfs"

IN_VCF="$VCFDIR/batch1234CN_DmagnaLRV01_merged1_10.vcf.gz"
FILT_VCF="$VCFDIR/batch1234CN_filtered.vcf.gz"

# ── Settings ───────────────────────────────────────────────────────
MIN_QUAL=30             # mimimum site quality
MIN_MAF=0.05            # minor allele frequency cutoff
MAX_MISSING=0.50        # maximum fraction of missing genotypes per site
MIN_DEPTH=5             # minimum mean sample depth
MIN_GQ=20               # minimum genotype quality
THREADS=8               # CPU threads

# -------------------------------------------------------------------------------
# Filter with bcftools 
#    site-level:
#	- keeps only biallelic SNPs (m2, M2, V indels)
#	- apply QUAL filter
#    genotype-level:
#	- mask genotypes that have low read depth and SNPs that have low GQ
#    site-level:
#	- remove sites based on % absence.
# -------------------------------------------------------------------------------

echo "================="

# Index the input VCF
echo "indexing the input vcf..."
tabix -p vcf "$IN_VCF"

# Filter 
echo "filtering..."
bcftools view \
    --threads $THREADS \
    -m2 \
    -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL" \
    -Ou \
    "$IN_VCF" \
| \
bcftools +setGT -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DEPTH) || (ALT!=\".\" && FMT/GQ<$MIN_GQ)" \
| \
bcftools view \
    --threads $THREADS \
    -i "F_MISSING<=$MAX_MISSING && MAF>=$MIN_MAF" \
    -Oz \
    -o "$FILT_VCF"

# Index the filtered VCF
echo "indexing the filtered vcf..."
bcftools index -t "$FILT_VCF"

echo "filtering is done"
