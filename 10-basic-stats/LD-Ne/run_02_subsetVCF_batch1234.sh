#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=VCF_subsetting
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-20 	# 20 populations

 
# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
cd /lustre1/scratch/363/vsc36396/LD_Ne/currentNe2

LIST_DIR="/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/09_basicstats/LD_Ne/populations_sample_lists"
VCF_DIR="vcfs"
mkdir -p "$VCF_DIR"

pops=(A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050)
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}
pop_samples="$LIST_DIR/${pop}_samples.txt"

IN="batch1234.filtered.vcf.gz"
OUT="$VCF_DIR/${pop}.filtered.vcf.gz"

	# Subset the VCF per population & index
echo "subsetting $IN ..."

bcftools view -S "$pop_samples" \
  -o "$OUT" -O z \
  --threads 10 "$IN"
bcftools index -t "$OUT"

echo "subsetting is done"

