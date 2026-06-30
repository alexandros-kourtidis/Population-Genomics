#!/usr/bin/env bash
# =============================================================================
# LD pruning & PCA with PLINK
# (PLINK2 uses a different file format and is stuck in 'beta' version)
# =============================================================================
#SBATCH --cluster=mindwell
#SBATCH --job-name=ld_prune
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/02_ld_prune.%j.log
#SBATCH -A lp_svbelleghem

set -euo pipefail
mkdir -p logs

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load PLINK/1.9
module load HTSlib/1.22.1-GCC-14.3.0

# ── Directories & Files ─────────────────────────────────────────────────────────
cd "/lustre1/scratch/363/vsc36396/treemix"
VCFDIR="vcfs"
OUTDIR="PLINK"

mkdir -p "$OUTDIR"

IN_VCF="$VCFDIR/batch1234CN_filtered.vcf.gz"
OUT="$OUTDIR/batch1234CN"

# ── Settings ───────────────────────────────────────────────────────
LD_WINDOW_SIZE=50                                      # Window size in SNPs
LD_STEP=10                                             # Step size in SNPs
LD_R2=0.1                                              # r² threshold for LD pruning
THREADS=8                                              # CPU threads

# perform linkage pruning - i.e. identify prune sites
echo "performing ld prunning..."
plink --vcf "$IN_VCF" --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise $LD_WINDOW_SIZE $LD_STEP $LD_R2 --threads $THREADS --out "$OUT"

# perform the pca and create .bed files for downstream analyses
echo "doing the pca..."
plink --vcf "$IN_VCF" --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract "${OUT}.prune.in"  --pca 100 --allow-no-sex  --threads $THREADS --out "$OUT" --make-bed  

# export pruned SNPs back to VCF for downstream use (e.g. TreeMix)
echo "exporting pruned VCF..."
plink --bfile "$OUT" \
      --allow-extra-chr \
      --recode vcf \
      --threads $THREADS \
      --out "${OUT}_ldpruned"

bgzip "${OUT}_ldpruned.vcf"
tabix -p vcf "${OUT}_ldpruned.vcf.gz"

echo "prunning and pca is done"

