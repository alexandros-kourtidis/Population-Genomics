#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/wice/batch
module load BCFtools/1.15.1-GCC-11.3.0
module load PLINK/1.9
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396

	# set up environmental variables
VCF=batch1234_DmagnaLRV01_merged1_10
OUT=batch1234_plink

	# Filtering 
bcftools view -i 'QUAL>=30 && MIN(GQ)>=30 && MIN(FMT/DP)>=10' -m2 -M2 --types snps ${VCF}.vcf.gz | \
bcftools filter -e 'F_MISSING>0.75' | \
bgzip --threads 10 > ${VCF}.filtered.vcf.gz
echo "Filtering by bcftools is done"

	# Index the filtered file
bcftools index -t ${VCF}.filtered.vcf.gz
echo "Indexing by bcftools is done"

	# perform linkage pruning & filter minimum allele frequency 
plink --vcf ${VCF}.filtered.vcf.gz --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --maf 0.02 --out $OUT

	# create bed files etc
plink --vcf ${VCF}.filtered.vcf.gz --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--extract ${OUT}.prune.in --make-bed --out $OUT  

echo "Plink has run"