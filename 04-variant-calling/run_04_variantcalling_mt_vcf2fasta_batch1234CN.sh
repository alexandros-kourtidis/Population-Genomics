#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=genotyping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 # Make sure the array size matches your sample count (D. magna: 10 chromosomes, 10+ scaffolds)
 
# Load the programs we will use
module load SAMtools/1.16.1-GCC-11.3.0
module load cluster/genius/amd
module load BCFtools/1.15.1-GCC-11.3.0
module load HTSlib/1.17-GCC-12.2.0
 
echo "================="
 
cd /lustre1/scratch/363/vsc36396/alignments_all/
IN_DIR=.
OUT_DIR="VCFs"
mkdir -p $OUT_DIR

REF=/lustre1/scratch/363/vsc36396/reference_LRV_01/DmagnaLRV01.fasta

	#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

	# Chromosome names and samples IDs
chrom=(Mt)
samples=($(cat /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/04_genotyping/batch1234CN_list_sequence.txt))
 
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
bcftools mpileup -O z --threads 16 -f $REF -r ${chrom[$ID]} $(echo $command) -a AD,DP | \
bcftools call -m -O z --threads 16 -f GQ,GP -o $OUT_DIR/CN_DmagnaLRV01_${chrom[$ID]}.vcf.gz

	# Index the VCF
bcftools index $OUT_DIR/CN_DmagnaLRV01_${chrom[$ID]}.vcf.gz

	# Activate conda for the vcf2fasta environment
source /vsc-hard-mounts/leuven-data/363/vsc36396/miniconda3/etc/profile.d/conda.sh
conda activate vcf2fasta

	# Extract a fasta file from the resulting VCF
vcf2fasta -f $REF $OUT_DIR/CN_DmagnaLRV01_${chrom[$ID]}.vcf.gz -o ${chrom[$ID]}_batch1234CN.fasta
 
echo "genotyping done"