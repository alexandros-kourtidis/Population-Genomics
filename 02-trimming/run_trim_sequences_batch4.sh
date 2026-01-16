#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-97 # Make sure the array size matches your sample count

module load Trimmomatic/0.39-Java-1.8.0_192
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(A5_GenM_C35 B4_OHZ_C48 C2_MO_C12 C2_MO_C18 C2_MO_C20 C2_MO_C22 C2_MO_C24 C2_MO_C26 C2_MO_C28 C2_MO_C32 C2_MO_C34 C2_MO_C36 C2_MO_C38 C2_MO_C42 C2_MO_C44 C4_BW_48630_C04 C4_BW_48630_C40 D4_BW_62256_C41 A2_BuSN_C08 A2_BuSN_C15 A2_BuSN_C16 A2_BuSN_C18 A2_BuSN_C20 A4_PL15_YEL_C14 A4_PL15_YEL_C16 A4_PL15_YEL_C18 A4_PL15_YEL_C19 A4_PL15_YEL_C22 A4_PL15_YEL_C28 A4_PL15_YEL_C30 A4_PL15_YEL_C32 A4_PL15_YEL_C36 A4_PL15_YEL_C38 A4_PL15_YEL_C44 A4_PL15_YEL_C46 A4_PL15_YEL_C48 A4_PL15_YEL_C50 B5_DA2_C01 B5_DA2_C02 B5_DA2_C06 B5_DA2_C23 B5_DA2_C26 B5_DA2_C27 B5_DA2_C29 B5_DA2_C30 C3_Ter1_C16 C3_Ter1_C18 C3_Ter1_C20 C3_Ter1_C22 C3_Ter1_C24 C3_Ter1_C28 C3_Ter1_C30 C3_Ter1_C32 C3_Ter1_C36 C3_Ter1_C38 C3_Ter1_C40 C3_Ter1_C45 C3_Ter1_C50 A2_BuSN_C22 A2_BuSN_C24 A2_BuSN_C26 A2_BuSN_C28 A2_BuSN_C32 A2_BuSN_C34 A2_BuSN_C36 A2_BuSN_C40 A2_BuSN_C44 A2_BuSN_C46 A2_BuSN_C50 A5_GenM_C06 A5_GenM_C14 A5_GenM_C18 A5_GenM_C22 A5_GenM_C50 B1_BKN1_C44 B1_BKN1_C50 C5_BW_36962_C01 C5_BW_36962_C26 C1_BlfN_C02 C1_BlfN_C14 C1_BlfN_C16 C1_BlfN_C18 C1_BlfN_C22 C1_BlfN_C24 C1_BlfN_C26 C1_BlfN_C28 C1_BlfN_C30 C1_BlfN_C34 C1_BlfN_C36 C1_BlfN_C38 C1_BlfN_C40 C1_BlfN_C44 C1_BlfN_C47 C3_Ter1_C42 C5_BW_36962_C32)

cd /lustre1/scratch/363/vsc36396/daphnia_reseq/batch4
 
#run trimmomatic for trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_unpaired.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_BGIadapters.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50
