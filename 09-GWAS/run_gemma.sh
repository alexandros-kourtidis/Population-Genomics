#!/bin/bash -l

#SBATCH --clusters=genius
#SBATCH --job-name=gemma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem

############################
# User-configurable inputs #
############################
PHENO_FILE="phenotype.txt"      # Header must include: FID, IID, Phenotype
PHENO_COL="Phenotype"           # Column header of the phenotype values in PHENO_FILE
IN="batch1234_plink"

#####################
# Environment setup #
#####################
module load tabix/0.2.6-GCCcore-6.4.0
module load PLINK/2.0.0-a.6.9-gfbf-2024a

# GEMMA environment
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh
conda activate gemma

#########################
# 0) Prep output dirs   #
#########################
cd /lustre1/scratch/363/vsc36396/gemma
mkdir -p kinship_matrix gemma_results logs

#########################################################
# 1) Clean the names in PLINK (A1_AnOudin_C02 -> A1C02) #
#########################################################
#plink \
#  --bfile $IN \
#  --allow-extra-chr \
#  --double-id \
#  --allow-no-sex \
#  --set-missing-var-ids @:# \
#  --update-ids new_headers.txt \
#  --make-bed \
#  --out ${IN}_renamed

##############################################################
# 2) Build a keep list of samples with NON-missing phenotype #
##############################################################
#   - Reads FID/IID/PHENO_COL from Phenotype.txt
#   - Skips NA/blank/-9 phenotypes
awk -v phcol="${PHENO_COL}" 'BEGIN{FS=OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="FID") f=i
    if($i=="IID") id=i
    if($i==phcol) p=i
  }
  if(!f || !id || !p){
    print "ERROR: Phenotype file must have FID, IID, and " phcol " headers." > "/dev/stderr"; exit 1
  }
  next
}
{
  val=$p
  gsub(/^[ \t]+|[ \t]+$/, "", val)
  if(val!="" && val!="NA" && val!="-9") print $f, $id
}' "${PHENO_FILE}" > keep_phenotyped.txt

if [ ! -s keep_phenotyped.txt ]; then
  echo "ERROR: No phenotyped samples found in ${PHENO_FILE} (column ${PHENO_COL})." >&2
  exit 1
fi

#######################################################
# 3) Subset to phenotyped samples, then add phenotype #
#######################################################
# Subset to the intersection of genotypes and valid phenotypes
plink \
  --bfile ${IN}_renamed \
  --keep keep_phenotyped.txt \
  --make-bed \
  --allow-extra-chr \
  --allow-no-sex \
  --out ${IN}_pheno_subset \
  > logs/step3_plink_keep.log 2>&1

# Add phenotype values into the .fam (column 6)
plink \
  --bfile ${IN}_pheno_subset \
  --pheno "${PHENO_FILE}" \
  --pheno-name "${PHENO_COL}" \
  --allow-extra-chr \
  --allow-no-sex \
  --make-bed \
  --out ${IN}_gwas_input \
  > logs/step3b_plink_add_pheno.log 2>&1

################################
# 4) Missingness (diagnostics) #
################################
plink \
  --bfile ${IN}_gwas_input \
  --missing \
  --allow-extra-chr \
  --allow-no-sex \
  --out missing_check \
  > logs/step4_missingness.log 2>&1

#################################
# 5) Kinship on phenotyped set  #
#################################
gemma \
  -bfile ${IN}_gwas_input \
  -gk 1 \
  -outdir kinship_matrix \
  -o "${IN}_gwas_input" \
  > logs/step6_gemma_gk.log 2>&1

################################
# 6) GEMMA LMM (no covariates) #
################################
gemma \
  -bfile "${IN}_gwas_input" \
  -k "kinship_matrix/${IN}_gwas_input.cXX.txt" \
  -lmm 4 \
  -outdir gemma_results \
  -o "${IN}_gemma_lmm_results_nocov" \
  > logs/step7_gemma_lmm.log 2>&1

echo "Done. Results in gemma_results/, kinship in kinship_matrix/, logs in logs/."
