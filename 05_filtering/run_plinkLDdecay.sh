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

	# perform linkage pruning & filter minimum allele frequency 
plink --vcf ${VCF}.filtered.vcf.gz --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --maf 0.05 --out $OUT

	# create bed files 
plink --vcf ${VCF}.filtered.vcf.gz --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# \
--extract ${OUT}.prune.in --make-bed --out $OUT  


plink --vcf ${VCF}.filtered.vcf.gz --double-id --allow-extra-chr --allow-no-sex \
	--make-bed --out data

	# Estimate pairwise LD (r²) in sliding windows
plink \
  --bfile data \
  --r2 \
  --ld-window 99999 \
  --ld-window-kb 500 \
  --ld-window-r2 0 \
  --out ld_decay

echo "Plink has run"