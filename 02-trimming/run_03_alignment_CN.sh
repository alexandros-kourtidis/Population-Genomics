#!/bin/bash -l
 
#SBATCH --cluster=genius
#SBATCH --job-name=alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=10:00:00
#SBATCH -A lp_svbelleghem
#SBATCH --array=1 # Make sure the array size matches your sample count (1-97)
 
#module loading
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.9-foss-2018a
module load Java

#This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))
 
# Sample IDs
samples=(CN_W1_1)
 
cd /lustre1/scratch/363/vsc36396/diversity_panel
mkdir CN_alignments
OUT_DIR="CN_alignments"
 
# first don't forget to index the reference genome - this is done only once
# Run BWA mapping - mem is the best choice for our DNB data
bwa mem -t 10 -M ~/SCRATCH/reference_LRV_01/DmagnaLRV01.fasta \
        $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_1_trimmed.fq.gz  \
        $(echo "${samples[ID]}")/$(echo "${samples[ID]}")_2_trimmed.fq.gz  \
| samtools view -bS > $OUT_DIR/$(echo "${samples[ID]}")_aligned.bam
 
# Filter using samtools
samtools view -f 0x02 -q 20 -b $OUT_DIR/$(echo "${samples[ID]}")_aligned.bam > $OUT_DIR/$(echo "${samples[ID]}")_aligned.filtered.bam
 
# Sort using samtools
samtools sort $OUT_DIR/$(echo "${samples[ID]}")_aligned.filtered.bam -o $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.bam
 
# Remove PCR duplicates
java -jar /vsc-hard-mounts/leuven-data/363/vsc36396/programs/picard_old.jar MarkDuplicates \
-I $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.bam \
-O $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-M $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.dup_metrics.txt \
-ASSUME_SORTED true
 
# Index using samtools
samtools index $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.nd.bam

# Remove intermediate files
# rm $OUT_DIR/$(echo "${samples[ID]}")_aligned.bam : this time don't execute this line because we need them to see the quality of the alignments with flagstat
 rm $OUT_DIR/$(echo "${samples[ID]}")_aligned.filtered.bam
 rm $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.bam
 
# Check the alignment quality with flagstat
samtools flagstat $OUT_DIR/$(echo "${samples[ID]}")_aligned.bam > $OUT_DIR/$(echo "${samples[ID]}")_alignment.flagstat.txt
samtools flagstat $OUT_DIR/$(echo "${samples[ID]}")_aligned.sorted.nd.bam > $OUT_DIR/$(echo "${samples[ID]}")_alignment.sorted.nd.flagstat.txt
echo "done"
 
###so the final BAM file that goes for genotyping is batch4_alignments/$(echo "${samples[ID]}")_aligned.sorted.nd.bam