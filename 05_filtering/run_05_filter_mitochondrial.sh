#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_Mt
#SBATCH --output=logs/filter_Mts.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8    # bcftoools only uses cores for compressing/decompressing so more would be a waste
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# =================================================================================
# In general: 3-step filtering with bcftools 
# ---------------------------------------------------------------------------------
#  1. site-level (we keep none for the mitochondrial here):
#	  - keeps only invariant sites (-V indels,mnps,snps)
#	    -NOT- apply QUAL filter
#     - apply total-DP filter (only upper limit - INFO/DP)
#	  - apply MQ filter
#  2. genotype-level:
#	  - mask genotypes with low/high read depth (FORMAT/DP)
#	    -NOT- mask SNPs that have low GQ
#  3. site-level:
#	  - remove sites that high MISS.
#	    -NOT- remove sites with low MAF
#  =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/filtering

INDIR="/lustre1/scratch/XXX/vscXXXXX/variantcalling/VCFs"
VCFDIR="VCFs"
SUMDIR="summaries"
mkdir -p "$VCFDIR" "$SUMDIR"

base="batch1234_Mt"
IN_VCF="${INDIR}/${base}.vcf.gz"

VCF1="${VCFDIR}/${base}.raw.vcf.gz"
VCF2="${VCFDIR}/${base}.invsnps.noMitoMiss.vcf.gz"
VCF3="${VCFDIR}/${base}.invsnps.noMitoMiss_DP.vcf.gz"
VCF4="${VCFDIR}/${base}.invsnps.noMitoMiss_DP_MISS.vcf.gz"

sum1="${SUMDIR}/${base}.raw.summary.txt"
sum2="${SUMDIR}/${base}.invsnps.noMitoMiss.summary.txt"
sum3="${SUMDIR}/${base}.invsnps.noMitoMiss_DP.summary.txt"
sum4="${SUMDIR}/${base}.invsnps.noMitoMiss_DP_MISS.summary.txt"

# ── Settings ───────────────────────────────────────────────────────

#MIN_QUAL=40             # mimimum site quality, no because we need the invariant sites
#MIN_MQ=40               # minimum mapping quality, I checked manually they all have 59-60
#MAX_TOTALDP=20000       # maximum total site-level depth (INFO/DP) - mean per-sample depth (20) x n samples (400) x 2.5

MIN_DP=150                # minimum sample depth (FORMAT/DP) - based on the mitochondrial DP distribution
MAX_DP=280				# maximum sample depth (FORMAT/DP) - based on the mitochondrial DP distribution
#MIN_GQ=30               # minimum genotype quality, Phred score - with GQ=30: chance for incorrect call = 0.1%

MAX_MISS=0.10           # maximum fraction of missing genotypes per site
#MIN_MAF=0.02            # minor allele frequency cutoff - 20 populations, MAF=0.05 removes alleles fixed in one population

CLEAN_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/batch1234_rename.txt"
RM_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/exclude_samples" # file with each sample to be excluded in a separate row

# -------------------------------------------------------------------------------
#  VCF1: Clean samples' names
# -------------------------------------------------------------------------------
echo "cleaning sample names of $IN_VCF..."
bcftools reheader -s "$CLEAN_SAMPLES" -o "$VCF1" "$IN_VCF" 

echo "indexing the $VCF1..."
bcftools index -t "$VCF1"

# Creating stats of the cleaned vcf file
echo "creating stats txt file for $VCF1..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF1" > "$sum1"


# -------------------------------------------------------------------------------
# VCF2: Indel removal and contaminated/low quality samples (from nuclear) removal
# -------------------------------------------------------------------------------
echo "filtering for $VCF2 ..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -M2 \
    -V indels,mnps \
    -S ^"$RM_SAMPLES" \
    -Oz \
    -o "$VCF2" \
    "$VCF1" 

echo "indexing the filtered vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF2" > "$sum2"


# -------------------------------------------------------------------------------
# VCF3: Genotype-level Masking
# -------------------------------------------------------------------------------
echo "filtering for $VCF3 ..."

bcftools +setGT -Oz -o "$VCF3" "$VCF2" -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP | FMT/DP>$MAX_DP)" 

echo "indexing the filtered vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF3"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF3..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF3" > "$sum3"


# -------------------------------------------------------------------------------
# VCF4: Site-level Missingness (no MAF) Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF4 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "F_MISSING<=$MAX_MISS" \
    -Oz \
    -o "$VCF4" \
    "$VCF3"

echo "indexing the filtered vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF4"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF4..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF4" > "$sum4"

echo "done"