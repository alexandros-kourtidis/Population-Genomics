#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=filter_invs
#SBATCH --output=logs/filter_invs.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8    # bcftoools only uses cores for compressing/decompressing so more would be a waste
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# =================================================================================
# In general: 3-step filtering with bcftools 
# ---------------------------------------------------------------------------------
#  1. site-level:
#	  - keeps only invariant sites(-M2, -V snps,indels,mnps)
#	  - apply QUAL filter
#     - apply total-DP filter (only upper limit - INFO/DP)
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
cd /lustre1/scratch/363/vsc36396/SVs

IN_DIR="/lustre1/scratch/363/vsc36396/filtering/VCFs"
VCFDIR="vcfs"
SUMDIR="summaries"
mkdir -p "$VCFDIR" "$SUMDIR"

base="batch1234_merged1_10"
IN_VCF="${IN_DIR}/${base}.raw.vcf.gz"
VCF1="${VCFDIR}/${base}.inv.QUALDPMQ.vcf.gz"
VCF2="${VCFDIR}/${base}.inv.QUALDPMQ_DP.vcf.gz"
VCF3="${VCFDIR}/${base}.inv.QUALDPMQ_DP_MISS.vcf.gz"

sum0="${SUMDIR}/${base}.raw.summary.txt"
sum1="${SUMDIR}/${base}.inv.QUALDPMQ.summary.txt"
sum2="${SUMDIR}/${base}.inv.QUALDPMQ_DP.summary.txt"
sum3="${SUMDIR}/${base}.inv.QUALDPMQ_DP_MISS.summary.txt"


# ── Settings ───────────────────────────────────────────────────────

MIN_QUAL=40             # mimimum site quality, Phred score - with 40, chance to be wrongly called =0.01%
MIN_MQ=40               # minimum mapping quality, Phred score - with 40, chance to be wrongly mapped =0.01%
MAX_TOTALDP=20000       # maximum total site-level depth (INFO/DP) - mean per-sample depth (20) x n samples (400) x 2.5

MIN_DP=8                # minimum sample depth (FORMAT/DP) - with DP=8, chance that we have seen both alleles = 99.2% with allele frequency=0.5
MAX_DP=80				# maximum sample depth (FORMAT/DP) - mean per-sample depth(20) x 4
#MIN_GQ=30               # minimum genotype quality, Phred score - with GQ=30: chance for incorrect call = 0.1%

MAX_MISS=0.10           # maximum fraction of missing genotypes per site
#MIN_MAF=0.02            # minor allele frequency cutoff - 20 populations, MAF=0.05 removes alleles fixed in one population

# Index the input vcf file
#bcftools index -t "$IN_VCF"
# Creating stats of the input vcf file
#echo "creating stats txt file for $IN_VCF ..."
#bcftools stats -s - "$IN_VCF" > "$sum0"
# Creating miss files
#echo "creating imiss files for $IN_VCF ..."
#vcftools --gzvcf "$IN_VCF" --missing-site --out "$miss0"
#vcftools --gzvcf "$IN_VCF" --missing-indv --out "$miss0"

# -------------------------------------------------------------------------------
# VCF1: Site-level Quality Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF1 ..."

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -M2 \
    -V indels,mnps,snps \
    -i "QUAL>=$MIN_QUAL && \
        MQ>=$MIN_MQ && \
        INFO/DP<=$MAX_TOTALDP" \
    -Oz \
    -o "$VCF1" \
    "$IN_VCF" 

echo "indexing the filtered vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF1"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF1..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF1" > "$sum1"


# -------------------------------------------------------------------------------
# VCF2: Genotype-level Masking
# -------------------------------------------------------------------------------
echo "filtering for $VCF2 ..."

bcftools +setGT -Oz -o "$VCF2" "$VCF1" -- \
	-t q \
	-n . \
    -i "(FMT/DP<$MIN_DP | FMT/DP>$MAX_DP)" 

echo "indexing the filtered vcf..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF2" > "$sum2"

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
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF3"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF3..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF3" > "$sum3"

echo "done"