#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=admixture
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=21-25

	# Load the programs
conda activate admixture
 
echo "================="
cd /lustre1/scratch/363/vsc36396/alignments_all/VCFs

	# File filtered for: callinng-genotype quality, depth, % in the population, 
	# indel removal, removing non-variant sites, linkage prunned, maf<0.05
FILE=batch1234_nmaf

	# Make a directory in the home directory
mkdir -p /lustre1/scratch/363/vsc36396/ADMIXTURE
cd /lustre1/scratch/363/vsc36396/ADMIXTURE

	# Generate the input file in plink format. Will use the one I did it for the pca
#plink --vcf /home/data/vcf/$FILE.vcf.gz --make-bed --out $FILE --allow-extra-chr

	# ADMIXTURE does not accept chromosome names that are not human chromosomes. Will just exchange the first column by 0.
#awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
#mv $FILE.bim.tmp $FILE.bim
	
	# Run ADMIXTURE for a series number of K
K=$SLURM_ARRAY_TASK_ID
admixture --cv $FILE.bed $K > log${K}.out

echo "ADMIXTURE is done"
