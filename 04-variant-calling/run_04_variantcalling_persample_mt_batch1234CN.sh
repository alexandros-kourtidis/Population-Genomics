#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --job-name=genotyping_single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-398  # Adjust this to match the number of samples (397 batch1234 + 1 CN_W1_1)

# Load required modules
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0

echo "================="

cd /lustre1/scratch/363/vsc36396/alignments_all/
IN_DIR=.
OUT_DIR="VCFs_individual"
mkdir -p $OUT_DIR

REF=/lustre1/scratch/363/vsc36396/reference_LRV_01/DmagnaLRV01.fasta
chrom="Mt"

# Get sample name from list
samples=($(cat /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/04_genotyping/batch1234CN_list_sequence.txt))
sample=${samples[$((SLURM_ARRAY_TASK_ID - 1))]}

# Run mpileup and call for single sample
bcftools mpileup -O z -f $REF -r $chrom ${sample}_aligned.sorted.nd.bam -a AD,DP | \
bcftools call -m -O z -f GQ,GP -o $OUT_DIR/${sample}_${chrom}.vcf.gz

# Index the VCF
bcftools index $OUT_DIR/${sample}_${chrom}.vcf.gz

echo "Single-sample VCF for $sample done"