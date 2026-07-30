#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=concat_invsnps
#SBATCH --output=logs/concat_invsnps.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8    # bcftoools only uses cores for compressing/decompressing so more would be a waste
#SBATCH --mem-per-cpu=2g
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# =================================================================================
# Concatenate the two vcfs to create a invariant-snps only vcf after separate filtering 
# =================================================================================
set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/363/vsc36396/filtering

VCFDIR="VCFs"
SUMDIR="summaries"
TMP_SORT="${VCFDIR}/tmp_sort_${SLURM_JOB_ID}"

mkdir -p "$VCFDIR" "$SUMDIR" "$TMP_SORT"

base="batch1234_merged1_10"
IN_VCF1="${VCFDIR}/${base}.inv.DPMQ_DP_MISS_noMitoMiss.vcf.gz"
IN_VCF2="${VCFDIR}/${base}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"

VCF1="${VCFDIR}/${base}.invsnp.filtered_sorted.vcf.gz"

sum1="${SUMDIR}/${base}.invsnp.filtered_sorted.summary.txt"

# ── Settings ───────────────────────────────────────────────────────

# -------------------------------------------------------------------------------
# VCF1: Site-level Quality Filter
# -------------------------------------------------------------------------------
echo "concatenating $IN_VCF1 and $IN_VCF2 ..."

bcftools concat \
    --threads ${SLURM_CPUS_PER_TASK} \
    -a \
    "$IN_VCF1" \
    "$IN_VCF2" \
    -Ou \
| \
bcftools sort \
    -T "$TMP_SORT" \
    -m 8G \
    -Oz \
    -o "$VCF1" \
    -
rm -rf "$TMP_SORT"

echo "indexing the concatenated $VCF1 ..."
bcftools index -t --threads ${SLURM_CPUS_PER_TASK} "$VCF1"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF1..."
bcftools stats -s - --threads ${SLURM_CPUS_PER_TASK} "$VCF1" > "$sum1"

echo "done"