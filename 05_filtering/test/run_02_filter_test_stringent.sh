#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_test
#SBATCH --output=logs/filtertest.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# =================================================================================
#  Checking the filters. For this we will only work with VARIANT SITES to reduce 
#   the working size and because they are the informative ones in most cases. 
# =================================================================================
# In general: 3-step filtering with bcftools 
# ---------------------------------------------------------------------------------
#  1. site-level:
#	  - keeps only biallelic SNPs (m2, M2, V indels)
#	  - apply QUAL filter
#     - apply total-DP filter (only upper limit)
#	  - apply MQ filter
#  2. genotype-level:
#	  - mask genotypes with low/high read depth (FORMAT/DP)
#	  - mask SNPs that have low GQ
#  3. site-level:
#	  - remove sites that high MISS.
#	  - remove sites with low MAF
#  =================================================================================
#   This filter test is to check the very hard filtering happening in the MISS step
#  because of the FORMAT/DP masking 
#  =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
module load VCFtools/0.1.16-GCC-10.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/filtering

IN_DIR="/lustre1/scratch/XXX/vscXXXXX/variantcalling/VCFs"
VCFDIR="test/vcfs"
SUMDIR="test/summaries"
mkdir -p "$VCFDIR" "$SUMDIR"

base="batch1234_scaffold_10"
IN_VCF="${IN_DIR}/${base}.vcf.gz"
VCF2="${VCFDIR}/${base}.stringent.vcf.gz"

sum0="${SUMDIR}/${base}.raw.summary.txt"
sum2="${SUMDIR}/${base}.stringent.summary.txt"


# ── Settings ───────────────────────────────────────────────────────

MIN_QUAL2=30             # mimimum site quality, Phred score - with QUAL=20: chance = 1%
MIN_MQ2=40               # minimum mapping quality, Phred score
MAX_TOTALDP=20000       # maximum total site-level depth (INFO/DP) - mean per-sample depth (20) x n samples (400) x 2.5

MIN_DP=8                # minimum sample depth (FORMAT/DP) - with DP=5, chance that we have seen both alleles = 93.8% 
MAX_DP=80				# maximum sample depth (FORMAT/DP) - mean per-sample depth(20) x 4
MIN_GQ=30               # minimum genotype quality, Phred score - with GQ=30: chance for incorrect call = 0.1%

MAX_MISS=0.10           # maximum fraction of missing genotypes per site
MIN_MAF=0.02           # minor allele frequency cutoff

# Index the input vcf file
#bcftools index -t "$IN_VCF"
# Creating stats of the input vcf file
#echo "creating stats txt file for $IN_VCF ..."
#bcftools stats -s - "$IN_VCF" > "$sum0"

# -------------------------------------------------------------------------------
# VCF2: MISS 0.10
# -------------------------------------------------------------------------------

echo "filtering for $VCF2 ..."
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL2 && \
        MQ>=$MIN_MQ2 && \
        INFO/DP<=$MAX_TOTALDP" \
    -Ou \
    "$IN_VCF" \
| \
bcftools +setGT -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP | FMT/DP>$MAX_DP) | \
        (ALT!=\".\" && FMT/GQ<$MIN_GQ)" \
| \
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "F_MISSING<=$MAX_MISS && MAF>=$MIN_MAF" \
    -Oz \
    -o "$VCF2"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - "$VCF2" > "$sum2"

# Creating miss files
echo "creating imiss files for $VCF2 ..."
vcftools --gzvcf $VCF2 --missing-site --out "${VCF2}_miss"
vcftools --gzvcf $VCF2 --missing-indv --out "${VCF2}_miss"


echo "done"