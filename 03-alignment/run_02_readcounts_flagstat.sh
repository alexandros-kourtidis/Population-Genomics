#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name count_average
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=count_average.%j.out
#SBATCH -A lp_svbelleghem

module load SAMtools/1.9-foss-2018a

# Directory containing the BAM and flagstat files
bam_dir="/lustre1/scratch/363/vsc36396/daphnia_reseq/batch4_alignments"

# Output file
output_file="/lustre1/scratch/363/vsc36396/daphnia_reseq/batch4_alignments/readcount_flagstat_batch4.txt"

# Initialize the output file with headers
echo -e "Filename\tTotal Reads" > "$output_file"

# Iterate over each flagstat file in the directory
for flagstat_file in "$bam_dir"/*.flagstat.txt; do
  # Extract the filename
  filename=$(basename "$flagstat_file" .flagstat.txt)

  # Parse the total reads from the flagstat file
  total_reads=$(awk 'NR==1 {print $1}' "$flagstat_file")

  # Write the result to the output file
  echo -e "$filename\t$total_reads" >> "$output_file"
done

echo "Total read counts have been written to $output_file."
