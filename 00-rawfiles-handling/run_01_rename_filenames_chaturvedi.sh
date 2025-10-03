#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=renaming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH -A lp_svbelleghem
 
 
# Directory containing the directories to rename
cd /lustre1/scratch/363/vsc36396/daphnia_reseq/DATA_Chaturvedi
 
#For making sure that the txt file created on Windows is used in a Unix-like environment:
dos2unix /vsc-hard-mounts/leuven-data/363/vsc36396/scripts/01_renaming/chaturvedi_list.txt
 
#Path to the renaming file
RENAME_FILE=/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/01_renaming/chaturvedi_list.txt
 
 
# Rename directories based on the text file Rename directories based on the text file
while IFS=$'\t' read -r new_name identifier; do
        # Rename the directory
        # mv S-${identifier} ${new_name}
 
        # Navigate to the newly renamed directory
        # cd ${new_name}
 
        # Rename the files within each directory
        mv ${identifier}_1.fastq.gz ${new_name}_1.fq.gz
        mv ${identifier}_2.fastq.gz ${new_name}_2.fq.gz
 
        # Go back to the parent directory
        # cd ..
 
done < "$RENAME_FILE"
 
echo "Renaming of directories and files completed."
 