#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=4-5 # Make sure the array size matches your sample count

module load Trimmomatic/0.39-Java-1.8.0_192
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(D3_BKLE5_C08 D3_BKLE5_C09 A3_Mech_C12 A3_Mech_C14 A3_Mech_C18 A3_Mech_C20 A3_Mech_C24 A3_Mech_C26 A3_Mech_C28 A3_Mech_C34 A3_Mech_C36 B2_OM2_C44 B4_OHZ_C20 B4_OHZ_C21 B4_OHZ_C22 B4_OHZ_C24 B4_OHZ_C26 B4_OHZ_C29 B4_OHZ_C31 C4_BW_48630_C05 D2_LRV_C38 D2_LRV_C40 A3_Mech_C38 B2_OM2_C14 B2_OM2_C18 B4_OHZ_C23 B4_OHZ_C25 B4_OHZ_C30 B4_OHZ_C32 B4_OHZ_C33 B4_OHZ_C37 B4_OHZ_C38 B4_OHZ_C39 B4_OHZ_C40 B4_OHZ_C43 B4_OHZ_C46 B5_DA2_C04 B5_DA2_C05 B5_DA2_C07 B5_DA2_C08 B5_DA2_C13 B5_DA2_C17 A3_Mech_C40 A3_Mech_C41 A3_Mech_C46 A3_Mech_C48 A3_Mech_C49 A5_GenM_C02 A5_GenM_C10 A5_GenM_C20 A5_GenM_C23 A5_GenM_C24 A5_GenM_C30 A5_GenM_C36 A5_GenM_C38 A5_GenM_C39 A5_GenM_C41 A5_GenM_C43 A5_GenM_C46 A5_GenM_C48 A5_GenM_C49 B2_OM2_C22 B2_OM2_C28 B2_OM2_C31 B2_OM2_C32 B2_OM2_C34 B2_OM2_C36 B2_OM2_C38 B2_OM2_C40 D1_CBOO6_C12 D1_CBOO6_C16 D1_CBOO6_C18 D1_CBOO6_C20 D1_CBOO6_C22 D1_CBOO6_C24 D1_CBOO6_C26 D1_CBOO6_C28 D1_CBOO6_C30 D1_CBOO6_C32 D1_CBOO6_C34 D1_CBOO6_C36 D1_CBOO6_C42 D1_CBOO6_C44 D1_CBOO6_C46 D3_BKLE5_C11 D3_BKLE5_C14 D3_BKLE5_C15 D3_BKLE5_C16 D3_BKLE5_C17 D3_BKLE5_C18 D3_BKLE5_C19 D3_BKLE5_C22 D3_BKLE5_C23 D3_BKLE5_C24 D3_BKLE5_C25 D3_BKLE5_C26 D3_BKLE5_C28 B4_OHZ_C50 D4_BW_62256_C16 D4_BW_62256_C43)
 
cd /lustre1/scratch/363/vsc36396/daphnia_reseq/batch3
 
#run trimmomatic for trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_unpaired.fq.gz \
$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_BGIadapters.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50
