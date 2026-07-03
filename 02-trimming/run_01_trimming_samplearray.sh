#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH -A lp_svbelleghem

module load cluster/genius/amd
module load Trimmomatic/0.39-Java-1.8.0_192

cd /lustre1/scratch/363/vsc36396/sf2/CN_W1_1

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
CN_W1_1_1.fq.gz CN_W1_1_2.fq.gz \
CN_W1_1_1_trimmed.fq.gz CN_W1_1_1_unpaired.fq.gz \
CN_W1_1_2_trimmed.fq.gz CN_W1_1_2_unpaired.fq.gz \
ILLUMINACLIP:/vsc-hard-mounts/leuven-data/363/vsc36396/scripts/02_trimming/TruSeqadapters_fastqc.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50

echo "trimming is done"
