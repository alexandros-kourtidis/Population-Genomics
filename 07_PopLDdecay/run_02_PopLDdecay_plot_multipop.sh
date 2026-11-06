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
plotscript_multipop="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/gitscripts/bin/Plot_MultiPop.pl"

# ── Directories & Files ─────────────────────────────────────────────────────────
PLOT_DIR1="plots/combined500kb_part1"
PLOT_DIR2="plots/combined500kb_part2"
mkdir -p  "${PLOT_DIR1}" "${PLOT_DIR2}"

filelist1="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/outfile_list_part1.txt"
filelist2="/user/leuven/XXX/vscXXXXX/DATA/scripts/07_PopLDdecay/outfile_list_part2.txt"

plotfile1="${PLOT_DIR1}/combined_part1.LDdecay500kb"
plotfile2="${PLOT_DIR2}/combined_part2.LDdecay500kb"

# -------------------------------------------------------------------------------
# Plot 1: LD decay for first batch
# -------------------------------------------------------------------------------
echo "Plotting LD decay for multiple populations part 1..."

perl "$plotscript_multipop" \
	-inList "${filelist1}" \
	-output "${plotfile1}" \
	-maxX 500
echo "PopLDdecay has run"

# -------------------------------------------------------------------------------
# Plot 2: LD decay for second batch
# -------------------------------------------------------------------------------
echo "Plotting LD decay for multiple populations part 2..."

perl "$plotscript_multipop" \
	-inList "${filelist2}" \
	-output "${plotfile2}" \
	-maxX 500
echo "PopLDdecay has run"