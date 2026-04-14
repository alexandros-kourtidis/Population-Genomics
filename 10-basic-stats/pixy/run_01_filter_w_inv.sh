#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=pixy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
 
# -------------------------------------------------------------------------------
#  Filtering VCF for pixy keeping non-variant sites: 
#   indel removal; callinng quality (QUAL, all sites); biallelic only; missingness in population;
#   genotype quality (GQ, only on variant sites); sequencing depth (DP, all sites); 
# -------------------------------------------------------------------------------

# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.18-GCC-12.3.0
module load SciPy-bundle/2023.07-iimkl-2023a

echo "================="

cd /lustre1/scratch/363/vsc36396/fst

	# Output directories
VCF_DIR="."
OUT_DIR="filtered_vcfs"
mkdir -p "$VCF_DIR" "$OUT_DIR" 

	# Filenames
IN_VCF="$VCF_DIR/batch1234_DmagnaLRV01_merged1_10.vcf.gz"
FILT_VCF="$OUT_DIR/batch1234.filtered_w_inv.vcf.gz"
GT_VCF="$OUT_DIR/batch1234_w_inv.filtered.gtonly.vcf.gz"

# -------------------------------------------------------------------------------
# 1) Filter with bcftools (reduce size but keep non-variant sites)
#   site-level:
#    - biallelic (<=2 alleles)
#    - removes indels, keep SNPs AND invariant sites with -V indels (!)
#    - apply QUAL filter
#   genotype-level:
#    - mask genotypes that have low read depth and SNPs that have low GQ
#   site-level:
#    - remove sites that are absent in 50% of samples or more.
# -------------------------------------------------------------------------------

echo "Step 1: bcftools filtering"

bcftools view --threads 10\
	-M2 \
	-V indels \
	-i 'QUAL>=10' \
	-Ou \
	"$IN_VCF" \
| \
bcftools +setGT --threads 10 -- \
	-t q \
	-n . \
	-i '(FMT/DP<8) || (ALT!="." && FMT/GQ<30)' \
| \
bcftools view --threads 10 \
	-i 'F_MISSING<=0.50' \
	-Oz \
	-o "$FILT_VCF"

bcftools index -t --threads 10 "$FILT_VCF"


# -------------------------------------------------------------------------------
# 2) Keep only GT in FORMAT
#    (drop all other FORMAT fields to polarise with polarizeVCFbyOutgroup.py)
# -------------------------------------------------------------------------------

echo "Step 2: retain only GT in FORMAT"

bcftools annotate --threads 10 \
	-x ^FORMAT/GT \
	-Oz -o "$GT_VCF" \
	"$FILT_VCF"
bcftools index -t --threads 10 "$GT_VCF"


# -------------------------------------------------------------------------------
# 3) Clean sample names in the final calls
# -------------------------------------------------------------------------------

#echo "Step 5: cleaning sample names"

#sed -i 's/_aligned\.sorted\.nd\.bam//g' "$CALLS_OUT"

#echo "Calls file created: $CALLS_OUT"
#echo "${pop} Done."
