#!/bin/bash -l
 
# -------------------------------------------------------------------------------
#   remove unnecessary samples 
# -------------------------------------------------------------------------------

#SBATCH --cluster=wice
#SBATCH --job-name=cleaning_VCF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --output=logs/01_clean.%j.log
#SBATCH -A lp_svbelleghem

# Environment setup
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# Variables
cd /user/leuven/XXX/vscXXXXX/SCRATCH/filtering
IN_VCF="batch1234chat_DmagnaLRV01_merged1_10.vcf.gz"
RM_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/chaturvedi_list.txt"
VCF1="batch1234.vcf.gz"

# -------------------------------------------------------------------------------
# Remove Chaturvedi's samples
# -------------------------------------------------------------------------------
echo "removing samples..."

bcftools view --threads 5 -S ^"$RM_SAMPLES" \
  -Oz -o "$VCF1" \
  "$IN_VCF"
bcftools index -t "$VCF1"

echo "done"