#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 # Make sure the array size matches your sample count (1-97)

module load cluster/genius/amd
module load Trimmomatic/0.39-Java-1.8.0_192
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(CN_W1_1)

cd /lustre1/scratch/363/vsc36396/diversity_panel
 
#run trimmomatic for trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_unpaired.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_trimming/TruSeqadapters_fastqc.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50

echo "trimming done!"