#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=01s
#SBATCH --output=logs/filter_invs.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1    # bcftoools query doesn't multithread
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem

# ── Environment ───────────────────────────────────────────────────────────
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0

# ── Scaffold ───────────────────────────────────────────────────────────
SCAFFOLD="scaffold_4"

# ── Directories & Files ─────────────────────────────────────────────────────────
cd /lustre1/scratch/XXX/vscXXXXX/SVs

VCFDIR="vcfs"
OUTDIR="01presences"
mkdir -p "$VCFDIR" "$OUTDIR"

base="batch1234_merged1_10"
VCF1="${VCFDIR}/${base}.inv.QUALDPMQ.vcf.gz"

out1="${OUTDIR}/${base}.inv.QUALDPMQ.${SCAFFOLD}.01.txt"

# ── Samples names ───────────────────────────────────────────────────────────
samples=$(bcftools query -l "$VCF1")

# get all invariant sites as 0 or 1 for a specific scaffold 
echo "extracting locations of invariant sites on $SCAFFOLD from $VCF1 ..."

{
  # header line
  printf "CHROM\tPOS\t%s\n" "$(echo "$samples" | tr '\n' '\t' | sed 's/\t$//')"

  # data lines
  bcftools query -r "$SCAFFOLD" \
  -f '%CHROM\t%POS[\t%GT]\n' "$VCF1" \
| awk 'BEGIN{OFS="\t"} {
    out = $1 OFS $2
    for (i = 3; i <= NF; i++) {
      out = out OFS ($i != "./." ? 1 : 0)
    }
    print out
  }'
} > "$out1"

echo "done"