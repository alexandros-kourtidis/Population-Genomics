#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=count_average
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=count_average.%j.out
#SBATCH -A lp_svbelleghem
 
# Output file
output="read_stats.tsv"
echo -e "File\t\t\tRead_Count\tAvg_Read_Length" > "$output"
 
# Loop through each directory within the directory you are using the command in
for dir in */; do
    sample=$(basename "$dir")  # Get directory name as sample name
 
    # Process each FASTQ file inside the directory
    for fastq in "$dir"/*{1,2,trimmed}.fq.gz; do
        # Get read count (total lines / 4)
        read_count=$(zcat "$fastq" | wc -l | awk '{print $1/4}')
 
        # Get average read length
        avg_length=$(zcat "$fastq" | awk '{if(NR%4==2) {count++; bases += length($0)} } END{if(count>0) print bases/count; else print "NA"}')
 
        # Write results to file
        echo -e "$(basename "$fastq")\t\t$read_count\t$avg_length" >> "$output"
    done
done
