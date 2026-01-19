#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=genotyping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-14 # Make sure the array size matches your sample count (D. magna: 10 chromosomes, 10+ scaffolds)
 
# Load the programs we will use
module load SAMtools/1.16.1-GCC-11.3.0
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
	# Store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

	# Paths and variables
cd /lustre1/scratch/363/vsc36396/alignments_all
mkdir -p /lustre1/scratch/363/vsc36396/alignments_all/VCFs

OUT_DIR=/lustre1/scratch/363/vsc36396/alignments_all/VCFs
INPUT="batch1234chat_DmagnaLRV01"
REF=/lustre1/scratch/363/vsc36396/reference_chaturvedi/reference.fasta

	# Chromosome names
chrom=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14)

	# Sample IDs
samples=($(cat /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/batch1234chat_list_sequence.txt))
 
	# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[@]}
do
        ALL_LIST="$ALL_LIST ${FILE}_aligned.sorted.nd.bam"
done
 
	# Use eval to ensure the full list is constructed correctly
eval command="\$ALL_LIST"

	# run mpileup
	# Non-variant sites are included.
	# Will include both allelic and total depths (AD and DP), and the output format will be GT:GQ:GP (genotype: genotype quality: genotype posterior probability)
bcftools mpileup -Ou --threads 16 -f $REF -r ${chrom[$ID]} $(echo $command) -a AD,DP | bcftools call -m -O z -f GQ,GP -o ${OUT_FOLDER}/${INPUT}_${chrom[$ID]}.vcf.gz
 
echo "genotyping done"
