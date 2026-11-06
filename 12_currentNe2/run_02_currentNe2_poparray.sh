#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=currentNe2
#SBATCH --output=logs/currentNe2.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=3		#20 populations

# -------------------------------------------------------------------------------
#  Running currentNe2 for calculation of the current Ne based on LD decay
#  Sources:
#	https://github.com/esrud/currentNe2
#  The required genetic genome size (in Morgans) calculation is from:
#	https://doi.org/10.1186/s12863-016-0445-7 (Dukic et al., 2016)
# -------------------------------------------------------------------------------
#  Environment variables
#	program path:
# 	/vsc-hard-mounts/leuven-data/XXX/vscXXXXX/programs/currentNe2/currentNe2/currentne2
# -------------------------------------------------------------------------------


# ── Environment ───────────────────────────────────────────────────────────
module load GCC
module load OpenMPI

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/currentNe2

IN_DIR="VCFs"
OUT_DIR="currentNe2_out/no_substructure"
mkdir -p "${OUT_DIR}"

populations=(A1_AnOudin A2_BuSN A3_Mech A4_PL15YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW48630 C5_BW36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW62256 D5_BW22050)
ID=$((SLURM_ARRAY_TASK_ID - 1))
pop="${populations[$ID]}"

VCF="${IN_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss_popMISSMAF.vcf.gz"
VCF_UNZIPPED="${IN_DIR}/${pop}.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss_popMISSMAF.vcf"
OUT="${OUT_DIR}/${pop}"

# ── Settings ───────────────────────────────────────────────────────
R_GENOME_SIZE=16
THREADS=${SLURM_CPUS_PER_TASK}

# -------------------------------------------------------------------------------
# 1. Unzip the vcf file
# -------------------------------------------------------------------------------
echo "Decompressing VCF file ..."
gunzip -c "$VCF" > "$VCF_UNZIPPED"

# -------------------------------------------------------------------------------
# 2. Run currentNe2 
# 	-x accounts for structure inside the populations
#	The last argument is genome size in Morgans
# -------------------------------------------------------------------------------
echo "Running currentNe2 for ${pop}..."
/vsc-hard-mounts/leuven-data/XXX/vscXXXXX/programs/currentNe2/currentNe2/currentne2 \
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
