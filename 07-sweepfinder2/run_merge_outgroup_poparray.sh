#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=merging
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-20		# 20 populations
 
# -------------------------------------------------------------------------------
#  Merging each population vcf with the outgroup file, using bcftools merge
#   in a population-level array 
# -------------------------------------------------------------------------------

# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.18-GCC-12.3.0

echo "================="

cd /lustre1/scratch/363/vsc36396/sf2

	# Population IDs
pops=( \
	A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM \
	B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 \
	C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 \
	D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050\
)

	# Get the population ID for this task
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}

	# Directories
VCF_DIR="population_VCFs"
mkdir -p "$VCF_DIR"

	# Filenames
OUTGROUP_VCF="CN_VCFs/CN_DmagnaLRV01_merged1_10.vcf.gz"
IN_VCF="$VCF_DIR/${pop}.vcf.gz"
OUT_VCF="${pop}_CN.vcf.gz"
if [[ ! -f "$IN_VCF" ]]; then
  echo "ERROR: Input VCF not found: $IN_VCF" >&2
  exit 1
fi

	# Index your VCF files (if not already done so)
tabix -p vcf ${OUTGROUP_VCF}
tabix -p vcf ${IN_VCF}

	# Merge the VCFs of the each population with the outgroup sample.
bcftools merge --threads 10 ${OUTGROUP_VCF} ${IN_VCF} -o ${VCF_DIR}/${OUT_VCF} -O z

echo "${pop} Done."