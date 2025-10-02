#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-5 # Make sure the array size matches your sample count

module load Trimmomatic/0.39-Java-1.8.0_192
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(A1_AnOudin_C02_combined A1_AnOudin_C06_combined C5_BW_36962_C18_combined D4_BW_62256_C09_combined D4_BW_62256_C10_combined)

cd /lustre1/scratch/363/vsc36396/daphnia_reseq/duplicates
 
#run trimmomatic for trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_unpaired.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_BGIadapters.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50
