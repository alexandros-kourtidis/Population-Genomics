#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=vcf_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/genius/amd
module load VCFtools
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/sf2/VCFs
VCF="CN_DmagnaLRV01"

	# Merge all scaffold VCFs into a single VCF file.
bcftools concat --threads 10 -o ${VCF}_merged1_10.vcf.gz -O z ${VCF}_scaffold_1.vcf.gz ${VCF}_scaffold_2.vcf.gz ${VCF}_scaffold_3.vcf.gz ${VCF}_scaffold_4.vcf.gz ${VCF}_scaffold_5.vcf.gz ${VCF}_scaffold_6.vcf.gz ${VCF}_scaffold_7.vcf.gz ${VCF}_scaffold_8.vcf.gz ${VCF}_scaffold_9.vcf.gz ${VCF}_scaffold_10.vcf.gz

	#Calculate missing data for each individual/sample and prints it to an out.imiss file
vcftools --gzvcf ${VCF}_merged1_10.vcf.gz --missing-indv --out ${VCF}_merged1_10

echo "merging and missingess count is done!"