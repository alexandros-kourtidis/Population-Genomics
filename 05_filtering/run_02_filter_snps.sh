#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_test
#SBATCH --output=logs/filtertest.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8    # bcftoools only uses cores for compressing/decompressing so more would be a waste
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# =================================================================================
# In general: 3-step filtering with bcftools 
# ---------------------------------------------------------------------------------
#  1. site-level:
#	  - keeps only biallelic SNPs (-m2, -M2, -v snps)
#	  - apply QUAL filter
#     - apply total-DP filter (only upper limit - INFO/DP)
#	  - apply MQ filter
#  2. genotype-level:
#	  - mask genotypes with low/high read depth (FORMAT/DP)
#	  - mask SNPs that have low GQ
#  3. site-level:
#	  - remove sites that high MISS.
#	  - remove sites with low MAF
# 
# Additionally: Removal of samples with
#    - High missingness (>0.2)                                   --> A1C16, C2C06, C5C21
#    - Double mitochondrial haplotypes (possible contamination)  --> A1C14, D2C32
#  =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/filtering

VCFDIR="VCFs"
SUMDIR="summaries"
MISSDIR="missfiles"
mkdir -p "$VCFDIR" "$SUMDIR" "$MISSDIR"

base="batch1234_merged1_10"
IN_VCF="${VCFDIR}/${base}.raw.vcf.gz"
VCF1="${VCFDIR}/${base}.snp.QUALDPMQ.vcf.gz"
VCF2="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ.vcf.gz"
VCF3="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS.vcf.gz"
VCF4="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF.vcf.gz"
VCF5="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"
VCF6="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz"

sum0="${SUMDIR}/${base}.raw.summary.txt"
sum1="${SUMDIR}/${base}.snp.QUALDPMQ.summary.txt"
sum2="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ.summary.txt"
sum3="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS.summary.txt"
sum4="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF.summary.txt"
sum5="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.summary.txt"
sum6="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.summary.txt"

# ── Settings ───────────────────────────────────────────────────────

MIN_QUAL=40             # mimimum site quality, Phred score - with 40, chance to be wrongly called =0.01%
MIN_MQ=40               # minimum mapping quality, Phred score - with 40, chance to be wrongly mapped =0.01%
MAX_TOTALDP=20000       # maximum total site-level depth (INFO/DP) - mean per-sample depth (20) x n samples (400) x 2.5

MIN_DP=8                # minimum sample depth (FORMAT/DP) - with DP=8, chance that we have seen both alleles = 99.2% with allele frequency=0.5
MAX_DP=80				# maximum sample depth (FORMAT/DP) - mean per-sample depth(20) x 4
MIN_GQ=30               # minimum genotype quality, Phred score - with GQ=30: chance for incorrect call = 0.1%

MAX_MISS=0.10           # maximum fraction of missing genotypes per site
MIN_MAF=0.02            # minor allele frequency cutoff - 20 populations, MAF=0.05 is too much since it removes alleles fixed in one population

RM_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/exclude_samples" # file with each sample to be excluded in a separate row

# Index the input vcf file
bcftools index -t "$IN_VCF"
# Creating stats of the input vcf file
echo "creating stats txt file for $IN_VCF ..."
bcftools stats -s - "$IN_VCF" > "$sum0"


# -------------------------------------------------------------------------------
# VCF1: Site-level Quality Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF1 ..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m2 -M2 \
    -v snps \
    -i "QUAL>=$MIN_QUAL && \
        MQ>=$MIN_MQ && \
        INFO/DP<=$MAX_TOTALDP" \
    -Oz \
    -o "$VCF1" \
    "$IN_VCF" 

echo "indexing the filtered vcf..."
bcftools index -t "$VCF1"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF1..."
bcftools stats -s - "$VCF1" > "$sum1"

# -------------------------------------------------------------------------------
# VCF2: Genotype-level Masking
# -------------------------------------------------------------------------------
echo "filtering for $VCF2 ..."

bcftools +setGT -Oz -o "$VCF2" "$VCF1" -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP | FMT/DP>$MAX_DP) | \
        (ALT!=\".\" && FMT/GQ<$MIN_GQ)" 

echo "indexing the filtered vcf..."
bcftools index -t "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - "$VCF2" > "$sum2"

# -------------------------------------------------------------------------------
# VCF3: Site-level Missingness (no MAF) Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF3 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "F_MISSING<=$MAX_MISS" \
    -Oz \
    -o "$VCF3" \
    "$VCF2"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF3"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF3..."
bcftools stats -s - "$VCF3" > "$sum3"

# -------------------------------------------------------------------------------
# VCF4: Site-level Missingness and MAF Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF4 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "MAF>=$MIN_MAF" \
    -Oz \
    -o "$VCF4" \
    "$VCF3"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF4"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF4..."
bcftools stats -s - "$VCF4" > "$sum4"

# -------------------------------------------------------------------------------
# VCF5: Sample-level Missingness and Mitochondrial contamination removal from MISS
# -------------------------------------------------------------------------------
echo "filtering for $VCF5 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S ^"$RM_SAMPLES" \
    -Oz \
    -o "$VCF5" \
    "$VCF3"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF5"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF5..."
bcftools stats -s - "$VCF5" > "$sum5"

# -------------------------------------------------------------------------------
# VCF6: Sample-level Missingness and Mitochondrial contamination removal from MISSMAFF
# -------------------------------------------------------------------------------
echo "filtering for $VCF6 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -S ^"$RM_SAMPLES" \
    -Oz \
    -o "$VCF6" \
    "$VCF4"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF6"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF6..."
bcftools stats -s - "$VCF6" > "$sum6"

echo "done"