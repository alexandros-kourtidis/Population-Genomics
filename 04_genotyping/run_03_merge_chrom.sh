#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=vcf_merge
#SBATCH --output=logs/vcfmerge.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
  
cd /lustre1/scratch/XXX/vscXXXXX/variantcalling/VCFs
VCF="batch1234"
OUT="${VCF}_merged1_10.vcf.gz"
OUT_CLEAN="${VCF}_merged1_10.clean.vcf.gz"
CLEAN_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/batch1234_rename.txt"

	# Merge all scaffold VCFs into a single VCF file.
echo "merging the chromosome VCFs..."
bcftools concat --threads 24 -o "$OUT" -O z ${VCF}_scaffold_1.vcf.gz ${VCF}_scaffold_2.vcf.gz ${VCF}_scaffold_3.vcf.gz ${VCF}_scaffold_4.vcf.gz ${VCF}_scaffold_5.vcf.gz ${VCF}_scaffold_6.vcf.gz ${VCF}_scaffold_7.vcf.gz ${VCF}_scaffold_8.vcf.gz ${VCF}_scaffold_9.vcf.gz ${VCF}_scaffold_10.vcf.gz ${VCF}_Mt.vcf.gz

echo "indexing the VCFs..."
bcftools index -t "$OUT"

	# Clean sample names
echo "cleaning sample names..."
bcftools reheader -s "$CLEAN_SAMPLES" -o "$OUT_CLEAN" "$OUT" 
bcftools index -t "$OUT_CLEAN"

echo "done"