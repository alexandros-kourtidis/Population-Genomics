#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=renaming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-36 # Make sure the array size matches your sample count (1-36)

module load cluster/wice/batch_sapphirerapids
module load FastQC
 
# Directory containing the directories to rename
cd /lustre1/scratch/363/vsc36396/daphnia_reseq/chaturvedi
 
#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(Chaturvedi_OHZbottom_B1 Chaturvedi_OHZbottom_B2 Chaturvedi_OHZbottom_B11 Chaturvedi_OHZmedium_M2 Chaturvedi_OHZmedium_M3 Chaturvedi_OHZmedium_M4 Chaturvedi_OHZmedium_M5 Chaturvedi_OHZmedium_M6 Chaturvedi_OHZmedium_M7 Chaturvedi_OHZmedium_M8 Chaturvedi_OHZbottom_B3 Chaturvedi_OHZmedium_M9 Chaturvedi_OHZmedium_M10 Chaturvedi_OHZmedium_M11 Chaturvedi_OHZmedium_M12 Chaturvedi_OHZtop_T1 Chaturvedi_OHZtop_T2 Chaturvedi_OHZtop_T3 Chaturvedi_OHZtop_T4 Chaturvedi_OHZtop_T5 Chaturvedi_OHZtop_T6 Chaturvedi_OHZbottom_B4 Chaturvedi_OHZtop_T7 Chaturvedi_OHZtop_T8 Chaturvedi_OHZtop_T9 Chaturvedi_OHZtop_T10 Chaturvedi_OHZtop_T11 Chaturvedi_OHZtop_T12 Chaturvedi_OHZbottom_B5 Chaturvedi_OHZbottom_B6 Chaturvedi_OHZbottom_B7 Chaturvedi_OHZbottom_B8 Chaturvedi_OHZbottom_B9 Chaturvedi_OHZbottom_B10)

fastqc $(echo "${samples[ID]}")_1.fq.gz
 
echo "fastqc done"
 
