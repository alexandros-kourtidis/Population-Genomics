#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=vcf_rename
#SBATCH --output=logs/vcfmerge.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
 
# Load the programs we will use
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
  
VCF_DIR="/lustre1/scratch/XXX/vscXXXXX/filtering/VCFs"
cd "$VCF_DIR" || exit 1

CLEAN_SAMPLES="/user/leuven/XXX/vscXXXXX/DATA/scripts/05_filtering/CN_rename.txt"

shopt -s nullglob
VCF_FILES=("$VCF_DIR"/*.vcf.gz)

if (( ${#VCF_FILES[@]} == 0 )); then
	echo "No VCF files found in $VCF_DIR" >&2
	exit 1
fi

# Clean sample names in every VCF and keep each original filename.
for VCF in "${VCF_FILES[@]}"; do
	TEMP_VCF="${VCF%.vcf.gz}.renaming.tmp.vcf.gz"

	echo "Cleaning sample names: $(basename "$VCF")"
	if bcftools reheader -s "$CLEAN_SAMPLES" -o "$TEMP_VCF" "$VCF"; then
		mv "$TEMP_VCF" "$VCF"
		bcftools index -t -f "$VCF"
	else
		rm -f "$TEMP_VCF"
		echo "Failed to rename $VCF" >&2
		exit 1
	fi
done

echo "done"