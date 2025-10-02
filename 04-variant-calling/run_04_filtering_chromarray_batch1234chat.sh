#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=genotyping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-12 # Make sure the array size matches your sample count (D. magna: 10 chromosomes, 10+ scaffolds)
 
# Load the programs we will use
module load BWA
module load SAMtools/1.16.1-GCC-11.3.0
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
module load VCFtools
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/alignments_all/VCFs

	#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

	# Chromosome names
chrom=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14)

	#Calculate missing data for each individual/sample and prints it to an out.imiss file
vcftools --gzvcf batch1234chat_DmagnaLRV01_${chrom[$ID]}.vcf.gz --missing-indv --out batch1234chat_DmagnaLRV01_${chrom[$ID]}.vcf.gz 

 
	# Filtering in 4 levels: 1) callinng-genotype quality, depth and % of the population [for most analyses]; 2) indel removal [for analyses that can't handle them]; 3) removing non-variant sites 4) thinning [for PCA]
	# Now masked so that I only filter the merged files (different script).

	#1. Non-variant removal
	#maybe not necessary since the next step does it automatically...
#vcftools --gzvcf /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.noindel.vcf.gz --recode --min-alleles 2 \
#--stdout | bgzip > /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.noindel.varonly.vcf.gz

	#2. Quality filter:
	#CAREFUL! vcftools can not recognise non-variant sites' quality scores, and this step will actually remove the non-variant sites [if you use the GQ filter]
	# e.g.: simon martin's Fst scans needs the quality filtering to be done within it, so that it retains the non-variant sites
#vcftools --gzvcf /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.vcf.gz --recode --minQ 30 --minGQ 30 --minDP 10 --max-missing 0.25 \
#--stdout | bgzip > /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.vcf.gz

	#3. Indel removal:
#vcftools --gzvcf /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.vcf.gz --recode --remove-indels \
#--stdout | bgzip > /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.noindel.vcf.gz

	#4. Thinning:
#vcftools --gzvcf /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.noindel.varonly.vcf.gz --recode --thin 500 \
#--stdout | bgzip > /lustre1/scratch/363/vsc36396/alignments_batch123Chat/VCFs/batch123Chat_${chrom[$ID]}.filtered.noindel.varonly.thinned.vcf.gz


echo "done"