#!/bin/bash -l
 
#SBATCH --cluster=wice
#SBATCH --job-name=genotyping
#SBATCH --output=logs/01_genotyping.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-2 # Make sure the array size matches your sample count 

# Load the programs we will use
module load cluster/wice/batch
module load BCFtools/1.18-GCC-12.3.0
module load HTSlib/1.17-GCC-12.2.0
 
	# Store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

	# Paths and variables
cd /lustre1/scratch/363/vsc36396/variantcalling
IN_DIR=alignments
OUT_DIR="VCFs"
mkdir -p $OUT_DIR
VCF="batch1234"
REF="ref_genome/DmagnaLRV01.fasta"

	# Chromosome names and sample IDs
chrom=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14 Mt)
samples=($(cat /user/leuven/363/vsc36396/DATA/scripts/04_genotyping/lists/batch1234_list_sequence.txt))
 
	# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[@]}
do
        ALL_LIST="$ALL_LIST ${IN_DIR}/${FILE}_aligned.sorted.nd.bam"
done
 
	# Use eval to ensure the full list is constructed correctly
eval command="\$ALL_LIST"

# -------------------------------------------------------------------------------
	# run mpileup
	# Non-variant sites are included.
	# Will include both allelic and total depths (AD and DP), and the output format will be GT:GQ:GP (genotype: genotype quality: genotype posterior probability)
# -------------------------------------------------------------------------------
echo "running mpileup and variant calling..."
bcftools mpileup -Ou --threads 16 -f $REF -r ${chrom[$ID]} $(echo $command) -a AD,DP | \
bcftools call -m -O z --threads 16 -f GQ,GP -o ${OUT_DIR}/${VCF}_${chrom[$ID]}.vcf.gz

# -------------------------------------------------------------------------------
	# Index the VCFs
# -------------------------------------------------------------------------------
echo "indexing the VCFs..."
bcftools index -t $OUT_DIR/${VCF}_${chrom[$ID]}.vcf.gz

echo "done"