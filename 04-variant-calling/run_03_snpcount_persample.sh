#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=snp_count
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/sf2/VCFs
VCF="batch1234CN_DmagnaLRV01_merged1_10.vcf.gz"
VCF_NAME=$(basename "$VCF" .vcf.gz)

# Loop over all samples and, for each one, count sites where that sample's GT is non-reference (not 0/0, 0, or missing).
# Achieved by subsetting to one sample, filtering GT, and counting remaining variant rows with wc -l.
for s in $(bcftools query -l "$VCF"); do
  n=$(bcftools view -s "$s" -Ou "$VCF" \
      | bcftools view -i 'GT!="0/0" && GT!="0|0" && GT!="0" && GT!="./." && GT!="."' -H \
      | wc -l)
  echo -e "${s}\t${n}"
done > per_sample_counts_${VCF_NAME}.tsv


echo "snp count is done!"