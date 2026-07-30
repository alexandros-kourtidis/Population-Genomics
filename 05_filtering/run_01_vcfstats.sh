#!/bin/bash -l
 
# -------------------------------------------------------------------------------
#  Cleaning VCF:
#   create stat files for the created vcfs per genotype and per sample
#   plot them
# -------------------------------------------------------------------------------

#SBATCH --cluster=wice
#SBATCH --job-name=cleaning_VCFvcfstats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --output=logs/01_vcfstats.%j.log
#SBATCH -A lp_svbelleghem

# Environment setup
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
module load matplotlib/3.10.5-iimkl-2025b

# Variables
cd /user/leuven/363/vsc36396/SCRATCH/filtering
VCF="batch1234.clean.vcf.gz"

# 0. View the header of the vcf
# bcftools view -h "$IN_VCF"

# -------------------------------------------------------------------------------
# 1. Create stats of the vcf file
# -------------------------------------------------------------------------------
echo "creating stats txt files..."
bcftools stats -s - "$VCF" > "${VCF}_summary.txt"

echo " - ${VCF}_summary.txt :"
grep "^SN" "${VCF}_summary.txt"

# -------------------------------------------------------------------------------
# 2. Plot the stats file
# -------------------------------------------------------------------------------echo "plotting the stats..."
plot-vcfstats -p plots -P batch1234.clean.vcf.gz_summary.txt

# Check heterozygoisity
grep "^PSC" grep batch1234.clean.vcf.gz_summary.txt | awk '{print $3, $6}' > het_per_sample.txt
sort -k2 -n het_per_sample.txt | tail -10 > het_per_sample_top10.txt

echo "done"