#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=currentNe2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem

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

IN_DIR="."
OUT_DIR="currentNe2_out"
mkdir -p ${OUT_DIR}

VCF="$IN_DIR/batch1234.filtered.vcf.gz"
VCF_UNZIPPED="$IN_DIR/batch1234.filtered.vcf"
OUT="$OUT_DIR/batch1234"

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
	-t 8 \
	-o $OUT \
	$VCF_UNZIPPED \
	16

echo "currentNe2 has run"
