#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=filtering_parsing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-20 	# 20 populations
 
# Load the programs we will use
module load cluster/genius/amd
module load VCFtools
module load HTSlib/1.17-GCC-12.2.0
module load SciPy-bundle/2023.07-iimkl-2023a
 
echo "================="
cd /lustre1/scratch/363/vsc36396/sweepfinder

	# Filtering per population for Sweepfinder2 keeping non-variant sites: 
	# indel removal; callinng quality (QUAL, all sites); biallelic only; missingness in population;
	# genotype quality (GQ, only on non-variant sites); sequencing depth (DP, all sites); 

	# Population IDs
pops=(A1_AnOudin A2_BuSN A3_Mech A4_PL15_YEL A5_GenM B1_BKN1 B2_OM2 B3_ZW B4_OHZ B5_DA2 C1_BlfN C2_MO C3_Ter1 C4_BW_48630 C5_BW_36962 D1_CBOO6 D2_LRV D3_BKLE5 D4_BW_62256 D5_BW_22050)

	# Get the population ID for this task
pop=${pops[$((SLURM_ARRAY_TASK_ID - 1))]}

	# Output directories
VCF_DIR="VCFs"
OUT_DIR="final_calls"
mkdir -p "$VCF_DIR"
mkdir -p "$OUT_DIR"

	# 1: Filter with vcftools
#vcftools --gzvcf $VCF_DIR/${pop}.vcf.gz \
#	--recode --stdout \
#	--remove-indels \
#	--minQ 10 \
#	--max-alleles 2 \
#	--max-missing 0.25 \
#	| bgzip > $VCF_DIR/${pop}.filtered.vcf.gz

#	echo "Filtered vcf files created: $VCF_DIR/${pop}.filtered.vcf.gz"

	# 2: Parse VCF with custom script
python parseVCF.py \
	--gtf flag=GQ min=30 gtTypes=Het \
	--gtf flag=GQ min=30 gtTypes=HomAlt \
	--gtf flag=DP min=10 \
	--skipIndels \
	-i $VCF_DIR/${pop}.filtered.vcf.gz \
	| python filterGenotypes.py -if phased -of diplo > $OUT_DIR/${pop}.filtered.diplo.calls

	# 3: Clean the samples' names
sed -i 's/_aligned\.sorted\.nd\.bam//g' $OUT_DIR/${pop}.filtered.diplo.calls	

	echo "Calls file created: $OUT_DIR/${pop}.filtered.diplo.calls"



	#Filtering for Sweepfinder with bcftools only- didn't work to include non-variant sites
#bcftools view -i 'QUAL>=30' batch1234_DmagnaLRV01_merged1_10.vcf.gz\
#  | bcftools view -e 'TYPE="indel"' \
#  | bcftools filter -e 'FMT/GQ<30 || FMT/DP<10' \
#  | bcftools view -M2 \
#  | bcftools view -i 'F_MISSING<=0.25' \
#  -o batch1234_DmagnaLRV01_merged1_10.sweepfinderfiltered.vcf.gz -O z

#echo "filtering is done"