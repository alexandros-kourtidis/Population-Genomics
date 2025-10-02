#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-105 # Make sure the array size matches your sample count

module load Trimmomatic/0.39-Java-1.8.0_192
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(A1_AnOudin_C02 A1_AnOudin_C04 A1_AnOudin_C06 A1_AnOudin_C09 A1_AnOudin_C10 A1_AnOudin_C12 A1_AnOudin_C14 A1_AnOudin_C16 A1_AnOudin_C18 A1_AnOudin_C20 A1_AnOudin_C22 A1_AnOudin_C24 A1_AnOudin_C26 A1_AnOudin_C28 A1_AnOudin_C30 A1_AnOudin_C32 A1_AnOudin_C35 A2_BuSN_C02 A2_BuSN_C04 A2_BuSN_C10 A2_BuSN_C12 A3_Mech_C02 A3_Mech_C04 A3_Mech_C06 A3_Mech_C08 A3_Mech_C10 A4_PL15_YEL_C02 A4_PL15_YEL_C04 A4_PL15_YEL_C07 A4_PL15_YEL_C08 A4_PL15_YEL_C10 A4_PL15_YEL_C12 B1_BKN1_C02 B1_BKN1_C04 B1_BKN1_C06 B1_BKN1_C08 B1_BKN1_C10 B2_OM2_C02 B2_OM2_C04 B2_OM2_C06 B2_OM2_C08 B2_OM2_C10 B3_ZW_C02 B3_ZW_C06 B3_ZW_C08 B3_ZW_C12 C1_BlfN_C04 C1_BlfN_C06 C1_BlfN_C08 C1_BlfN_C10 C1_BlfN_C12 C2_MO_C02 C2_MO_C04 C2_MO_C06 C2_MO_C08 C2_MO_C10 C3_Ter1_C02 C3_Ter1_C04 C3_Ter1_C06 C3_Ter1_C07 C3_Ter1_C10 C3_Ter1_C12 C5_BW_36962_C03 C5_BW_36962_C05 C5_BW_36962_C06 C5_BW_36962_C07 C5_BW_36962_C08 C5_BW_36962_C09 C5_BW_36962_C10 C5_BW_36962_C11 C5_BW_36962_C13 C5_BW_36962_C14 C5_BW_36962_C15 C5_BW_36962_C17 C5_BW_36962_C18 C5_BW_36962_C20 C5_BW_36962_C21 C5_BW_36962_C22 C5_BW_36962_C25 D1_CBOO6_C02 D1_CBOO6_C04 D1_CBOO6_C06 D1_CBOO6_C08 D1_CBOO6_C10 D2_LRV_C02 D2_LRV_C04 D2_LRV_C06 D2_LRV_C08 D2_LRV_C10 D4_BW_62256_C02 D4_BW_62256_C03 D4_BW_62256_C04 D4_BW_62256_C05 D4_BW_62256_C09 D4_BW_62256_C10 D4_BW_62256_C12 D4_BW_62256_C13 D4_BW_62256_C18 D4_BW_62256_C19 D4_BW_62256_C32 D4_BW_62256_C40 A1_AnOudin_C02_Nucleospin A1_AnOudin_C06_Nucleospin D4_BW_62256_C09_Nucleospin D4_BW_62256_C10_Nucleospin)
 
cd /lustre1/scratch/363/vsc36396/daphnia_reseq/batch1
 
#run trimmomatic for trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_unpaired.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_BGIadapters.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50
