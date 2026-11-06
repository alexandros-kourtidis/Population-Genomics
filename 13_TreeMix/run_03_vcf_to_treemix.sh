#!/usr/bin/env bash
# =============================================================================
# Convert an LD-pruned VCF to TreeMix format
# Intended to use the LD-pruned dataset produced by alex/02_LDprunning_PCA.sh
# =============================================================================
#SBATCH --cluster=wice
#SBATCH --job-name=treemix_convert
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/03_vcf_to_treemix.%j.log
#SBATCH -A lp_svbelleghem

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load SciPy-bundle/2025.07-iimkl-2025b

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/treemix/

IN_DIR="plink_out"
OUT_DIR="treemix_input"

mkdir -p "$OUT_DIR"

INPUT_VCF="${IN_DIR}/batch1234CN.snps.MISS_plink10kb.vcf.gz"
POPMAP="/user/leuven/XXX/vscXXXXX/DATA/scripts/13_TreeMix/population_map.txt"
OUT_FRQ="${OUT_DIR}/treemix.frq.gz"

SCRIPT="/user/leuven/XXX/vscXXXXX/DATA/scripts/13_TreeMix/run_03_vcf_to_treemix.py"
PYTHON_BIN="${PYTHON_BIN:-${PYTHON:-python3}}"

# -------------------------------------------------------------------------------
#  Run the vcf_to_treemix  
# -------------------------------------------------------------------------------
# sanity checks 
if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "Input VCF not found: ${INPUT_VCF}" >&2
    exit 1
fi
if [[ ! -f "${POPMAP}" ]]; then
    echo "Population map not found: ${POPMAP}" >&2
    exit 1
fi

echo "Input VCF: ${INPUT_VCF}"
echo "Population map: ${POPMAP}"
echo "TreeMix output: ${OUT_FRQ}"

echo "Creating output directories..."
mkdir -p "$(dirname "${OUT_FRQ}")"

#  convert the LD-pruned VCF to TreeMix format 
echo "Converting LD-pruned VCF to TreeMix format..."
"$PYTHON_BIN" "${SCRIPT}" \
    --vcf "${INPUT_VCF}" \
    --popmap "${POPMAP}" \
    --out "${OUT_FRQ}"

echo "TreeMix input written to ${OUT_FRQ}"
