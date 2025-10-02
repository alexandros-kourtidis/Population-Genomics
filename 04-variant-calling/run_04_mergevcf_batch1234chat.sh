#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=vcf_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/genius/amd
module load VCFtools
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/alignments_all/VCFs

	# Merge all scaffold VCFs into a single VCF file.
bcftools concat --threads 16 -o batch1234chat_DmagnaLRV01_merged1_10.vcf.gz -O z batch1234chat_DmagnaLRV01_scaffold_1.vcf.gz batch1234chat_DmagnaLRV01_scaffold_2.vcf.gz batch1234chat_DmagnaLRV01_scaffold_3.vcf.gz batch1234chat_DmagnaLRV01_scaffold_4.vcf.gz batch1234chat_DmagnaLRV01_scaffold_5.vcf.gz batch1234chat_DmagnaLRV01_scaffold_6.vcf.gz batch1234chat_DmagnaLRV01_scaffold_7.vcf.gz batch1234chat_DmagnaLRV01_scaffold_8.vcf.gz batch1234chat_DmagnaLRV01_scaffold_9.vcf.gz batch1234chat_DmagnaLRV01_scaffold_10.vcf.gz

	#Calculate missing data for each individual/sample and prints it to an out.imiss file
vcftools --gzvcf batch1234chat_DmagnaLRV01_merged1_10.vcf.gz --missing-indv --out batch1234chat_DmagnaLRV01_merged1_10.vcf.gz

echo "merging and missingess count is done!"