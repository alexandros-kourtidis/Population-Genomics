#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=genotyping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-14 # Make sure the array size matches your sample count (D. magna: 10 chromosomes, 10+ scaffolds)
 
# Load the programs we will use
module load SAMtools/1.16.1-GCC-11.3.0
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="

	# Directories
cd /lustre1/scratch/363/vsc36396/sf2
mkdir -p /lustre1/scratch/363/vsc36396/sf2/VCFs

	# Variables
IN_DIR="CN_alignments"
OUT_DIR="VCFs"
REF=/lustre1/scratch/363/vsc36396/lrv_ref/DmagnaLRV01.fasta

	#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

	# Chromosome names and samples IDs
chrom=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14)
samples=(CN_W1_1)
 
	# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[@]}
do
        ALL_LIST="$ALL_LIST $IN_DIR/${FILE}_aligned.sorted.nd.bam"
done
 
	# Use eval to ensure the full list is constructed correctly
eval command="\$ALL_LIST"

	# run mpileup
	# Non-variant sites are included.
	# Will include both allelic and total depths (AD and DP), and the output format will be GT:GQ:GP (genotype: genotype quality: genotype posterior probability)
bcftools mpileup -Ou --threads 10 -f $REF -r ${chrom[$ID]} $(echo $command) -a AD,DP | bcftools call --threads 10 -A -m -O z -f GQ,GP -o $OUT_DIR/CN_DmagnaLRV01_${chrom[$ID]}.vcf.gz
 
echo "genotyping done"