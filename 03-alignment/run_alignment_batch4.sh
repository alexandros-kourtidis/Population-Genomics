#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=2-97 # Make sure the array size matches your sample count (1-97)
 
#module loading
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.9-foss-2018a
module load Java

#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

cd /lustre1/scratch/363/vsc36396/daphnia_reseq
mkdir batch4_alignments

# Variables
samples=(A5_GenM_C35 B4_OHZ_C48 C2_MO_C12 C2_MO_C14 C2_MO_C16 C2_MO_C18 C2_MO_C20 C2_MO_C22 C2_MO_C24 C2_MO_C26 C2_MO_C28 C2_MO_C32 C2_MO_C34 C2_MO_C36 C2_MO_C38 C2_MO_C42 C2_MO_C44 C4_BW_48630_C04 C4_BW_48630_C40 D4_BW_62256_C41 A2_BuSN_C08 A2_BuSN_C15 A2_BuSN_C16 A2_BuSN_C18 A2_BuSN_C20 A4_PL15_YEL_C14 A4_PL15_YEL_C16 A4_PL15_YEL_C18 A4_PL15_YEL_C19 A4_PL15_YEL_C22 A4_PL15_YEL_C28 A4_PL15_YEL_C30 A4_PL15_YEL_C32 A4_PL15_YEL_C36 A4_PL15_YEL_C38 A4_PL15_YEL_C44 A4_PL15_YEL_C46 A4_PL15_YEL_C48 A4_PL15_YEL_C50 B5_DA2_C01 B5_DA2_C02 B5_DA2_C06 B5_DA2_C23 B5_DA2_C26 B5_DA2_C27 B5_DA2_C29 B5_DA2_C30 C3_Ter1_C16 C3_Ter1_C18 C3_Ter1_C20 C3_Ter1_C22 C3_Ter1_C24 C3_Ter1_C28 C3_Ter1_C30 C3_Ter1_C32 C3_Ter1_C36 C3_Ter1_C38 C3_Ter1_C40 C3_Ter1_C45 C3_Ter1_C50 A2_BuSN_C22 A2_BuSN_C24 A2_BuSN_C26 A2_BuSN_C28 A2_BuSN_C32 A2_BuSN_C34 A2_BuSN_C36 A2_BuSN_C40 A2_BuSN_C44 A2_BuSN_C46 A2_BuSN_C50 A5_GenM_C06 A5_GenM_C14 A5_GenM_C18 A5_GenM_C22 A5_GenM_C50 B1_BKN1_C44 B1_BKN1_C50 C5_BW_36962_C01 C5_BW_36962_C26 C1_BlfN_C02 C1_BlfN_C14 C1_BlfN_C16 C1_BlfN_C18 C1_BlfN_C22 C1_BlfN_C24 C1_BlfN_C26 C1_BlfN_C28 C1_BlfN_C30 C1_BlfN_C34 C1_BlfN_C36 C1_BlfN_C38 C1_BlfN_C40 C1_BlfN_C44 C1_BlfN_C47 C3_Ter1_C42 C5_BW_36962_C32)
REF="~/SCRATCH/reference_chaturvedi/reference.fasta"

#  don't forget to index the reference genome with bwa and samtools before!
# Run BWA mapping - mem is the best choice for our DNB data
bwa mem -t 10 -M $REF \
        /lustre1/scratch/363/vsc36396/daphnia_reseq/batch4/$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz  \
        /lustre1/scratch/363/vsc36396/daphnia_reseq/batch4/$(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz  \
| samtools view -bS > batch4_alignments/$(echo "${samples[ID]}")_aligned.bam
 
# Filter using samtools
samtools view -f 0x02 -q 20 -b batch4_alignments/$(echo "${samples[ID]}")_aligned.bam > batch4_alignments/$(echo "${samples[ID]}")_aligned.filtered.bam
 
# Sort using samtools
samtools sort -@ 10 batch4_alignments/$(echo "${samples[ID]}")_aligned.filtered.bam -o batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.bam

# Remove PCR duplicates
java -jar /vsc-hard-mounts/leuven-data/363/vsc36396/programs/picard_old.jar MarkDuplicates \
-I ~/SCRATCH/daphnia_reseq/batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.bam \
-O ~/SCRATCH/daphnia_reseq/batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-M ~/SCRATCH/daphnia_reseq/batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.dup_metrics.txt \
-ASSUME_SORTED true
 
# Index using samtools
samtools index batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.nd.bam

# Remove intermediate files
# rm batch4_alignments/$(echo "${samples[ID]}")_aligned.bam : this time don't execute this line because we need them to see the quality of the alignments with flagstat
 rm batch4_alignments/$(echo "${samples[ID]}")_aligned.filtered.bam
 rm batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.bam
 
# Check the alignment quality with flagstat
samtools flagstat batch4_alignments/$(echo "${samples[ID]}")_aligned.bam > batch4_alignments/$(echo "${samples[ID]}")_alignment.flagstat.txt
samtools flagstat batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.nd.bam > batch4_alignments/$(echo "${samples[ID]}")_alignment.sorted.nd.flagstat.txt
echo "done"
 
###so the final BAM file that goes for genotyping is batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.nd.bam
