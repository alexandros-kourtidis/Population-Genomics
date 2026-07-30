#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_test
#SBATCH --output=logs/filtertest.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
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

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/363/vsc36396/filtering

IN_DIR="/lustre1/scratch/363/vsc36396/variantcalling/VCFs"
VCFDIR="test/vcfs"
SUMDIR="test/summaries"
mkdir -p $VCFDIR $SUMDIR

base="batch1234_scaffold_10"
IN_VCF="${IN_DIR}/${base}.vcf.gz"
VCF1="${VCFDIR}/${base}.snpQUAL_filtered.vcf.gz"
VCF2="${VCFDIR}/${base}.snpQUALDPMQ_filtered.vcf.gz"
VCF3="${VCFDIR}/${base}.snpQUALDPMQ_DPGQ_filtered.vcf.gz"
VCF4="${VCFDIR}/${base}.snpQUALDPMQ_DPGQ_MISS_filtered.vcf.gz"

sum0="${SUMDIR}/${base}.raw.summary.txt"
sum1="${SUMDIR}/${base}.snpQUAL_filtered.summary.txt"
sum2="${SUMDIR}/${base}.snpQUALDPMQ_filtered.summary.txt"
sum3="${SUMDIR}/${base}.snpQUALDPMQ_DPGQ_filtered.summary.txt"
sum4="${SUMDIR}/${base}.snpQUALDPMQ_DPGQ_MISS_filtered.summary.txt"

# ── Settings ───────────────────────────────────────────────────────

MIN_QUAL=20             # mimimum site quality, Phred score - with QUAL=20: chance = 1%
MIN_MQ=30               # minimum mapping quality, Phred score
MAX_TOTALDP=32000       # maximum total site-level depth (INFO/DP) - mean per-sample depth (20) x n samples (400) x 4

MIN_DP=8                # minimum sample depth (FORMAT/DP) - with DP=5, chance that we have seen both alleles = 93.8% 
MAX_DP=80				# maximum sample depth (FORMAT/DP) - mean per-sample depth(20) x 4
MIN_GQ=30               # minimum genotype quality, Phred score - with GQ=30: chance for incorrect call = 0.1%

MAX_MISS=0.20           # maximum fraction of missing genotypes per site
MIN_MAF=0.05            # minor allele frequency cutoff

# Index the input vcf file
bcftools index -t "$IN_VCF"
# Creating stats of the input vcf file
echo "creating stats txt file for $IN_VCF ..."
bcftools stats -s - "$IN_VCF" > "$sum0"

# -------------------------------------------------------------------------------
# VCF1: Only QUAL filter 
# -------------------------------------------------------------------------------

echo "filtering for $VCF1 ..."
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 \
    -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL" \
    -Oz \
    -o "$VCF1" \
	"$IN_VCF"
echo "indexing the filtered vcf..."
bcftools index -t "$VCF1"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF1..."
bcftools stats -s - "$VCF1" > "$sum1"

# -------------------------------------------------------------------------------
# VCF2: QUAL+INFO/DP+MQ filter 
# -------------------------------------------------------------------------------

echo "filtering for $VCF2 ..."
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 \
    -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL && \
        MQ>=$MIN_MQ && \
        INFO/DP<=$MAX_TOTALDP"  \
    -Oz \
    -o "$VCF2" \
	"$IN_VCF"
echo "indexing the filtered vcf..."
bcftools index -t "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - "$VCF2" > "$sum2"

# -------------------------------------------------------------------------------
# VCF3: QUAL+INFO/DP+MQ & FORMAT/DP+GQ filter 
# -------------------------------------------------------------------------------

echo "filtering for $VCF3 ..."
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL && \
        MQ>=$MIN_MQ && \
        INFO/DP<=$MAX_TOTALDP" \
    -Ou \
    "$IN_VCF" \
| \
bcftools +setGT -Oz -o "$VCF3" -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP || FMT/DP>$MAX_DP) || \
        (ALT!=\".\" && FMT/GQ<$MIN_GQ)"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF3"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF3..."
bcftools stats -s - "$VCF3" > "$sum3"

# -------------------------------------------------------------------------------
# VCF4: QUAL+INFO/DP+MQ & FORMAT/DP+GQ & MISS filter 
# -------------------------------------------------------------------------------
echo "filtering for $VCF4 ..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL && \
        MQ>=$MIN_MQ && \
        INFO/DP<=$MAX_TOTALDP" \
    -Ou \
    "$IN_VCF" \
| \
bcftools +setGT -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP || FMT/DP>$MAX_DP) || \
        (ALT!=\".\" && FMT/GQ<$MIN_GQ)" \
| \
bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "F_MISSING<=$MAX_MISS && MAF>=$MIN_MAF" \
    -Oz \
    -o "$VCF4"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF4"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF4..."
bcftools stats -s - "$VCF4" > "$sum4"

echo "done"