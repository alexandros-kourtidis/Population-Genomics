#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/genius/amd
module load PLINK/1.9
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/PCA

# set up environmental variables
VCF=/lustre1/scratch/363/vsc36396/PCA/batch1234_DmagnaLRV01_merged1_10.PCAfiltered.vcf.gz

# perform linkage pruning - i.e. identify prune sites
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --maf 0.05 --out batch1234_nmaf

# perform the pca
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract batch1234_nmaf.prune.in  --pca 100 --allow-no-sex  --out batch1234_nmaf --make-bed  

echo "plink pca is done"