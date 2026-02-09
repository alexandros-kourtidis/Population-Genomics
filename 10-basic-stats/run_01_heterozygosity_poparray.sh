#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=heteroz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-20 	# 20 populations
 
# -------------------------------------------------------------------------------
#	Environment setup
# -------------------------------------------------------------------------------
module load cluster/genius/amd
module load VCFtools/0.1.16-GCC-10.3.0

echo "================="

# -------------------------------------------------------------------------------
#	Variables assignment
# -------------------------------------------------------------------------------

cd /lustre1/scratch/363/vsc36396/sf2/tmp_vcfs

	# Population 
pops=(A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050)

	# Get the population ID for this task
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}

	# Directories and files
OUT_DIR="het_stats"

IN_FILE="${pop}_CN.filtered.gtonly.vcf.gz"
OUT_FILE="$OUT_DIR/${pop}_CN.filtered.gtonly.het_per_sample"

mkdir -p "$OUT_DIR" 

# -------------------------------------------------------------------------------
#	Run heterozygosity per sampple scans 
# -------------------------------------------------------------------------------
#	O(HOM): Observed homozygous genotypes
#	E(HOM): Expected homozygous genotypes
#	N_SITES:Sites genotyped for this sample
#	F:	Inbreeding coefficient
#		>0 = inbred/excess homozygotes; 
#		<0 = excess heterozygotes
# -------------------------------------------------------------------------------

echo "Running heterozygosity test for ${pop}..."

time vcftools --gzvcf "$IN_FILE" --het --out "$OUT_FILE"

echo "Done"

