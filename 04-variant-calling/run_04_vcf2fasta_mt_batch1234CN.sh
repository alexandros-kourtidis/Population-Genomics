#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=vcf2fasta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
  
echo "================="
 
cd /lustre1/scratch/363/vsc36396/VCFs
VCF="CN_DmagnaLRV01_Mt"

REF=/lustre1/scratch/363/vsc36396/reference_LRV_01/DmagnaLRV01.fasta
GFF=/lustre1/scratch/363/vsc36396/reference_LRV_01/Daphnia_magna_LRV0_1_genes_with_Mt.gff3

	# Activate conda for the vcf2fasta environment
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh
conda activate vcf2fasta

	# Extract a fasta file from the resulting VCF
vcf2fasta -f $REF -v ${VCF}.vcf.gz -o ${VCF}.fasta
 
echo "vcf2fasta is done"