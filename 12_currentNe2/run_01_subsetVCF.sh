#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=VCF_subsetting
#SBATCH --output=logs/VCF_subsetting.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-20 	# 20 populations

 
# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
 
# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/currentNe2

LIST_DIR="/user/leuven/XXX/vscXXXXX/DATA/scripts/12_currentNe2/populations_sample_lists"
IN_DIR="/lustre1/scratch/XXX/vscXXXXX/VCFs"
OUT_DIR="VCFs"
mkdir -p "$OUT_DIR"

pops=(A1_AnOudin A2_BuSN A3_Mech A4_PL15YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW48630 C5_BW36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW62256 D5_BW22050)
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}
pop_samples="$LIST_DIR/${pop}_samples.txt"

IN="${IN_DIR}/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"
VCF1="${OUT_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"
VCF2="${OUT_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss_popMISSMAF.vcf.gz"

sum1="${OUT_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.summary.txt"
sum2="${OUT_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss_popMISSMAF.summary.txt"

# ── Settings ───────────────────────────────────────────────────────

MAX_MISS=0.10           # maximum fraction of missing genotypes per site
MIN_MAF=0.02            # minor allele frequency cutoff - 20 populations, MAF=0.05 is too much since it removes alleles fixed in one population

# -------------------------------------------------------------------------------
#  1. Subsetting 
# -------------------------------------------------------------------------------
#echo "subsetting ${IN} for ${pop}..."
#bcftools view -S "$pop_samples" \
#  --threads ${SLURM_CPUS_PER_TASK} \
#  -o "$VCF1" \
#  -O z \
#  "$IN"

#echo "indexing the filtered vcf..."
#bcftools index -t "$VCF1"

## Creating stats of the filtered vcf file
#echo "creating stats txt file for $OUT..."
#bcftools stats -s - "$OUT" > "$sum1"

# -------------------------------------------------------------------------------
# 2. Re-apply site-level Missingness and MAF Filter
# -------------------------------------------------------------------------------
echo "filtering for $VCF2 ..."     

bcftools view \
    --threads ${SLURM_CPUS_PER_TASK} \
    -i "F_MISSING<=$MAX_MISS && MAF>=$MIN_MAF" \
    -Oz \
    -o "$VCF2" \
    "$VCF1"

echo "indexing the filtered vcf..."
bcftools index -t "$VCF2"

# Creating stats of the filtered vcf file
echo "creating stats txt file for $VCF2..."
bcftools stats -s - "$VCF2" > "$sum2"

echo "done"
