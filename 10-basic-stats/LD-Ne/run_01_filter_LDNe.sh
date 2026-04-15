#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=pixy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
 
# -------------------------------------------------------------------------------
#  Filtering VCF for LD-based calculations of Ne: 
#   biallelic snps only; callinng quality (QUAL); missingness in population;
#   genotype quality (GQ); sequencing depth (DP); 
# -------------------------------------------------------------------------------

# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.18-GCC-12.3.0
module load SciPy-bundle/2023.07-iimkl-2023a

echo "================="

cd /lustre1/scratch/363/vsc36396/LD_Ne

	# Output directories
VCF_DIR="/lustre1/scratch/363/vsc36396/filtering_vcfs"
OUT_DIR="/lustre1/scratch/363/vsc36396/LD_Ne/currentNe2"
mkdir -p "$VCF_DIR" "$OUT_DIR" 

	# Filenames
IN_VCF="$VCF_DIR/batch1234_DmagnaLRV01_merged1_10.vcf.gz"
FILT_VCF="$OUT_DIR/batch1234.filtered.vcf.gz"

# -------------------------------------------------------------------------------
# 1) Filter with bcftools
#   site-level:
#    - biallelic 
#    - keep only SNPs 
#    - apply QUAL filter
#   genotype-level:
#    - mask genotypes that have low read depth and low GQ
#   site-level:
#    - keep sites that are present in at least 75% of samples.
# -------------------------------------------------------------------------------

echo "bcftools filtering..."

bcftools view --threads 15\
	-m2 -M2 \
	-v snps \
	-i 'QUAL>=20' \
	-Ou \
	"$IN_VCF" \
| \
bcftools +setGT --threads 15 -- \
	-t q \
	-n . \
	-i '(FMT/DP<8) || (FMT/GQ<10)' \
| \
bcftools view --threads 15 \
	-i 'F_MISSING<=0.25' \
	-Oz \
	-o "$FILT_VCF"

bcftools index -t --threads 15 "$FILT_VCF"

# -------------------------------------------------------------------------------
# 2) Counts in the filtered file:
#	- samples
#	- variants
#	- (non-missing) genotypes
# -------------------------------------------------------------------------------

echo "writting counts summary..."

nsamples=$(bcftools query -l "$FILT_VCF" | wc -l)
nvariants=$(bcftools index -n "$FILT_VCF")
ngenotypes=$(bcftools query -f '[%GT\t]\n' "$FILT_VCF" \
  | grep -o -E '[0-9]\|[0-9]|[0-9]/[0-9]' \
  | wc -l)

cat <<EOF > "$OUT"
VCF file:        $FILT_VCF
Samples:         $nsamples
Variants:        $nvariants
Non-missing GTs: $ngenotypes
EOF

echo "Filtering completed."

