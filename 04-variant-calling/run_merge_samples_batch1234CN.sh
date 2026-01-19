#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=vcf_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/genius/amd
module load VCFtools
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396
VCF_1="batch1234_DmagnaLRV01_merged1_10"
VCF_2="sf2/VCFs/CN_DmagnaLRV01_merged1_10"
VCF="batch1234CN_DmagnaLRV01_merged1_10"
OUT_DIR="sf2/VCFs/"

	# Index your VCF files (if not already done so)
tabix -p vcf ${VCF_1}.vcf.gz
tabix -p vcf ${VCF_2}.vcf.gz

	# Merge the VCFs of different samples ( into a single VCF file.
bcftools merge --threads 16 ${VCF_1}.vcf.gz ${VCF_2}.vcf.gz -o ${OUT_DIR}${VCF}.vcf.gz -O z

	#Calculate missing data for each individual/sample and prints it to an out.imiss file
vcftools --gzvcf ${OUT_DIR}${VCF}.vcf.gz --missing-indv --out $OUT_DIR${VCF}

echo "merging and missingess count is done!"
