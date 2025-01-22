#!/bin/bash
#SBATCH --job-name=filenames
#SBATCH --nodes=1
#SBATCH --time=0-00:10:00
#SBATCH --account=evolgen
#SBATCH --mem-per-cpu=10G

# Loop through the folders replicate_1 to replicate_100
for i in {1..100}; do
    folder="replicate_$i"
    output_file="rep_${i}_files.csv"

    # Check if the folder exists
    if [ -d "$folder" ]; then
        # List the files and save the output to a CSV file
        ls "$folder" > "$output_file"
        echo "Processed $folder -> $output_file"
    else
        echo "Folder $folder does not exist. Skipping."
    fi
done