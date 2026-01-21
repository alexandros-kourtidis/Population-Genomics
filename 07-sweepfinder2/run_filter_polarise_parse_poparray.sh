#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=filtering_parsing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 	# 20 populations
 
# -------------------------------------------------------------------------------
#  Filtering, polarising and parsing VCF per population for Sweepfinder2 keeping non-variant sites: 
#   indel removal; callinng quality (QUAL, all sites); biallelic only; missingness in population;
#   genotype quality (GQ, only on non-variant sites); sequencing depth (DP, all sites); 
# -------------------------------------------------------------------------------

# Load the programs we will use
module load cluster/genius/amd
module load BCFtools/1.18-GCC-12.3.0
module load SciPy-bundle/2023.07-iimkl-2023a

echo "================="

cd /lustre1/scratch/363/vsc36396/sf2

	# Population IDs
pops=( \
	A1_test A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM \
	B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 \
	C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 \
	D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050\
)

	# Get the population ID for this task
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}

	# Output directories
VCF_DIR="population_VCFs"
OUT_DIR="final_calls"
TMP_DIR="tmp_vcfs" 
mkdir -p "$VCF_DIR" "$OUT_DIR" "$TMP_DIR"

	# Filenames
IN_VCF="$VCF_DIR/${pop}_CN.vcf.gz"
if [[ ! -f "$IN_VCF" ]]; then
  echo "ERROR: Input VCF not found: $IN_VCF" >&2
  exit 1
fi
FILT_VCF="$TMP_DIR/${pop}_CN.filtered.vcf.gz"
GT_VCF="$TMP_DIR/${pop}_CN.filtered.gtonly.vcf.gz"
POL_VCF="$TMP_DIR/${pop}_CN.polarized.vcf.gz"
CALLS_OUT="$OUT_DIR/${pop}_CN.filtered.diplo.calls"
OUTGROUP="CN_alignments/CN_W1_1_aligned.sorted.nd.bam"

# -------------------------------------------------------------------------------
# 1) Filter with bcftools (reduce size but keep non-variant sites)
#    - remove indels, require biallelic (<=2 alleles)
#    - keep SNPs + REF sites
#    - apply QUAL and missingness filter to variant sites only
#      (reference-only sites pass regardless of QUAL)
# -------------------------------------------------------------------------------

echo "[${pop}] Step 1: bcftools filtering"

bcftools view \
	-M2 \
	--types snps,ref \
	-i 'QUAL>=10 && F_MISSING<=0.75' \
	-Oz -o "$FILT_VCF" \
	"$IN_VCF"
bcftools index -t "$FILT_VCF"


# -------------------------------------------------------------------------------
# 2) Keep only GT in FORMAT
#    (drop all other FORMAT fields to polarise with polarizeVCFbyOutgroup.py)
# -------------------------------------------------------------------------------

echo "[${pop}] Step 2: retain only GT in FORMAT"

bcftools annotate \
	-x ^FORMAT/GT \
	-Oz -o "$GT_VCF" \
	"$FILT_VCF"
bcftools index -t "$GT_VCF"


# -------------------------------------------------------------------------------
# 3) Polarize VCF with outgroup(s) with polarizeVCFbyOutgroup.py
#    -ind needs the exact sample name as in the header of the vcf
#    polarizeVCFbyOutgroup.py reads GT only;
#    output file cannot be indexed because it is not bgzipped (just gzipped) 
# -------------------------------------------------------------------------------

echo "[${pop}] Step 3: polarizing by outgroup"

python polarizeVCFbyOutgroup.py \
	-vcf "$GT_VCF" \
	-out "$POL_VCF" \
	-ind "${OUTGROUP}" \
	-add 

# bcftools index -t "$POL_VCF" 

# -------------------------------------------------------------------------------
# 4) Parse and filter genotypes with custom scripts, then convert to diploid
#    Notes:
#    - --gtf GQ filters will be applied (as configured) by your parser
#    - --skipIndels is harmless here since indels are already removed
# -------------------------------------------------------------------------------

echo "[${pop}] Step 4: parsing & genotype filtering"

python parseVCF.py \
	--gtf flag=GQ min=30 gtTypes=Het \
	--gtf flag=GQ min=30 gtTypes=HomAlt \
	--gtf flag=DP min=10 \
	--skipIndels \
	-i "$POL_VCF" \
	| python filterGenotypes.py -if phased -of diplo \
	> "$CALLS_OUT"


# -------------------------------------------------------------------------------
# 5) Clean sample names in the final calls
# -------------------------------------------------------------------------------

echo "[${pop}] Step 5: cleaning sample names"

sed -i 's/_aligned\.sorted\.nd\.bam//g' "$CALLS_OUT"

echo "Calls file created: $CALLS_OUT"
echo "${pop} Done."