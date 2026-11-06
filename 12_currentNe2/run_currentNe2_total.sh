#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=currentNe2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
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
# 	program path:
# 	/vsc-hard-mounts/leuven-data/XXX/vscXXXXX/programs/currentNe2/currentNe2/currentne2
# -------------------------------------------------------------------------------
cd /lustre1/scratch/XXX/vscXXXXX/LD_Ne/currentNe2

module load GCC
module load OpenMPI

IN_DIR="vcfs"
OUT_DIR="currentNe2_out"
mkdir -p "$OUT_DIR"

VCF="$IN_DIR/batch1234.filtered.vcf.gz"
VCF_UNZIPPED="$IN_DIR/batch1234.filtered.vcf"
OUT="$OUT_DIR/batch1234"

R_GENOME_SIZE=16
THREADS=${SLURM_CPUS_PER_TASK}

# -------------------------------------------------------------------------------
# 1. Unzip the vcf file
# -------------------------------------------------------------------------------

echo "Decompressing VCF file ..."

if [[ ! -f "$VCF_UNZIPPED" ]]; then
    echo "Unzipping $VCF ..."
    gunzip -c "$VCF" > "$VCF_UNZIPPED" || exit 1
else
    echo "Unzipped VCF already exists, skipping decompression."
fi

# -------------------------------------------------------------------------------
# 2. Run currentNe2 
#	Using a subsample (100000 out of ~ 700000SNPs) to be computationally feasible
#	The last argument is genome size in Morgans
# -------------------------------------------------------------------------------
echo "Running currentNe2 ..."

/vsc-hard-mounts/leuven-data/XXX/vscXXXXX/programs/currentNe2/currentNe2/currentne2 \
	-s 100000 \
	-t "$THREADS" \
	-o "$OUT" \
	$VCF_UNZIPPED \
	"$R_GENOME_SIZE"


# -------------------------------------------------------------------------------
# 3. Delete unzipped vcf file 
# -------------------------------------------------------------------------------
#echo "Deleting unzipped vcf file ..."

#rm -f "$VCF_UNZIPPED"

echo "currentNe2 has run"
