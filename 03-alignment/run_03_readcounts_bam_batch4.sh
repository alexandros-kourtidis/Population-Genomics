#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name count_average
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=count_average.%j.out
#SBATCH -A lp_svbelleghem

module load SAMtools/1.9-foss-2018a

# Directory containing the BAM files
bam_dir="/lustre1/scratch/363/vsc36396/daphnia_reseq/batch4_alignments"

# Output file
output_file="/lustre1/scratch/363/vsc36396/daphnia_reseq/batch4_alignments/readcount_bam_batch4.txt"

# Iterate over each BAM file in the directory
for bam_file in "$bam_dir"/*.nd.bam; do
  # Count the reads using samtools
  count=$(samtools view -c "$bam_file")
  # Write the result to the output file
  echo -e "$(basename "$bam_file")\t$count" >> "$output_file"
done

echo "Read counts have been written to $output_file."