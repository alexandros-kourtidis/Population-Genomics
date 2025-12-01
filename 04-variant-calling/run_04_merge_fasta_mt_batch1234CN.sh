#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --job-name=fasta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

module load libdeflate/1.20-GCCcore-13.3.0
module load BCFtools
module load SAMtools

echo "================="
cd /lustre1/scratch/363/vsc36396/alignments_all/VCFs_individual
REF=/lustre1/scratch/363/vsc36396/reference_LRV_01/DmagnaLRV01_Mt.fasta
OUT_DIR="FASTA_output"
mkdir -p $OUT_DIR

	# List of sample VCFs
samples=($(ls *.vcf.gz))

	# Loop through each sample VCF
for vcf in "${samples[@]}"; do
 sample_name=$(basename "$vcf" .vcf.gz)

	# Generate consensus FASTA for variant sites only
 bcftools consensus \
  -f $REF \
  -a '_' \
  --mark-del '-' \
  --mark-ins lc \
  "$vcf" | \
 tr -d '[a-z]' | \
 sed "1s/.*/>$sample_name/" \
 > "$OUT_DIR/${sample_name}_indel.fasta"

done

echo "All sample FASTAs generated."

	# Combine all FASTAs into one file
cat $OUT_DIR/*indel.fasta > $OUT_DIR/combined_samples_indel.fasta

echo "Combined FASTA created: $OUT_DIR/combined_samples_indel.fasta"