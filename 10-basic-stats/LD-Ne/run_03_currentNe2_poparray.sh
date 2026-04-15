#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=currentNe2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-20		#20 populations

# -------------------------------------------------------------------------------
#  Running currentNe2 for calculation of the current Ne based on LD decay
#  Sources:
#	https://github.com/esrud/currentNe2
#  The required genetic genome size (in Morgans) calculation is from:
#	https://doi.org/10.1186/s12863-016-0445-7 (Dukic et al., 2016)
# -------------------------------------------------------------------------------
#  Environment variables
# -------------------------------------------------------------------------------
cd /lustre1/scratch/363/vsc36396/LD_Ne/currentNe2

module load GCC
module load OpenMPI
# program path:
# /vsc-hard-mounts/leuven-data/363/vsc36396/programs/currentNe2/currentNe2/currentne2

populations=( A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM \
		B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 \
		C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 \
		D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050)
ID=$((SLURM_ARRAY_TASK_ID - 1))
pop="${populations[$ID]}"


IN_DIR="vcfs"
OUT_DIR="currentNe2_out"
mkdir -p "${OUT_DIR}"

VCF="$IN_DIR/${pop}.filtered.vcf.gz"
VCF_UNZIPPED="$IN_DIR/${pop}.filtered.vcf"
OUT="$OUT_DIR/${pop}"

R_GENOME_SIZE=16
THREADS=8

# -------------------------------------------------------------------------------
# 1. Unzip the vcf file
# -------------------------------------------------------------------------------
echo "Decompressing VCF file ..."

gunzip -c "$VCF" > "$VCF_UNZIPPED"

# -------------------------------------------------------------------------------
# 2. Run currentNe2 
#	The second argument is genome size in Morgans
# -------------------------------------------------------------------------------
echo "Running currentNe2 ..."

/vsc-hard-mounts/leuven-data/363/vsc36396/programs/currentNe2/currentNe2/currentne2 \
	-t "$THREADS" \
	-o "$OUT" \
	"$VCF_UNZIPPED" \
	"$R_GENOME_SIZE"

# -------------------------------------------------------------------------------
# 2. Delete unzipped vcf file 
# -------------------------------------------------------------------------------
echo "Deleting unzipped vcf file ..."

rm -f "$VCF_UNZIPPED"

echo "currentNe2 has run"
