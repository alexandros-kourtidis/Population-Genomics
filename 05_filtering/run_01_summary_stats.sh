#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=sum_stats
#SBATCH --output=logs/sum_stats.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/filtering

VCFDIR="VCFs"
SUMDIR="summaries"

base="batch1234_merged1_10"

VCF1="${VCFDIR}/${base}.inv.DPMQ_DP_MISS_noMitoMiss.vcf.gz"
VCF2="${VCFDIR}/${base}.invsnp.filtered_sorted.vcf.gz"
VCF3="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"
VCF4="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.vcf.gz"

sum1="${SUMDIR}/${base}.inv.DPMQ_DP_MISS_noMitoMiss.summary.txt"
sum2="${SUMDIR}/${base}.invsnp.filtered_sorted.summary.txt"
sum3="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.summary.txt"
sum4="${SUMDIR}/${base}.snp.QUALDPMQ_DPGQ_MISSMAF_noMitoMiss.summary.txt"

het1="${base}.inv.DPMQ_DP_MISS_noMitoMiss.hetpersample.txt"

# -------------------------------------------------------------------------------
#  All the summary stats + per sample
# -------------------------------------------------------------------------------

echo "creating stats txt file for $VCF1..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF1" > "$sum1"

echo "creating stats txt file for $VCF2..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF2" > "$sum2"

echo "creating stats txt file for $VCF3..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF3" > "$sum3"

echo "creating stats txt file for $VCF4..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF4" > "$sum4"

# -------------------------------------------------------------------------------
#  Check heterozygosity per sample
# -------------------------------------------------------------------------------
grep "^PSC" "$sum1" | awk '{print $3, $6}' > "$het1"

echo "done"