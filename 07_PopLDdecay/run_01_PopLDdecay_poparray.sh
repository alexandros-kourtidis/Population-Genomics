#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=popLDdecay
#SBATCH --output=logs/poplddecay_300.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=21 	#20 populations + 1 allpops
#SBATCH --mem=16G

# =================================================================================
#  Running PopLDdecay
# =================================================================================

# ── Environment ───────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/poplddecay
conda activate PopLDdecay
plotscript_onepop="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/gitscripts/bin/Plot_OnePop.pl"
plotscript_multipop="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/gitscripts/bin/Plot_MultiPop.pl"

# ── Directories & Files ─────────────────────────────────────────────────────────
VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/VCFs"
POP_DIR="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/population level"
OUT_DIR="ld_decay_out"
PLOT_DIR="plots"
mkdir -p "${OUT_DIR}" "${PLOT_DIR}"

VCF="${VCF_DIR}/batch1234_merged1_10.snp.QUALDPMQ_DPGQ_MISS_noMitoMiss.vcf.gz"
pops=(A1_AnOudin A2_BuSN A3_Mech A4_PL15YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW48630 C5_BW36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW62256 D5_BW22050 allpops)
pop="${pops[$((SLURM_ARRAY_TASK_ID - 1))]}"
poplist="${POP_DIR}/${pop}_samples.txt"
outfile="${OUT_DIR}/${pop}.LDdecay"
plotfile="${PLOT_DIR}/${pop}.LDdecay.png"

# ── Settings ───────────────────────────────────────────────────────
# only maximum distance we change the rest for now we use the defaults
DIST=500

# -------------------------------------------------------------------------------
# 1. Run PopLDdecay 
# -------------------------------------------------------------------------------
echo "Running PopLDdecay for population ${pop}..."

PopLDdecay \
    -InVCF "${VCF}" \
    -OutStat "${outfile}" \
    -SubPop "${poplist}" \
    -MaxDist $DIST 

# -------------------------------------------------------------------------------
# 2. Plot LD decay for the single population
# -------------------------------------------------------------------------------
echo "Plotting LD decay for population ${pop}..."

perl "$plotscript_onepop" \
	-inFile "${outfile}.stat.gz" \
	-output "${plotfile}"

echo "PopLDdecay has run"
