#!/bin/bash

# This function appends all the fasta files found in a directory (provided as positional argument 1) in a single file, and add the filename as the first field in the header.

# Verify if the first argument (path to find the fastas) was provided
if [ -z "$1" ]; then
	echo "Provide the path to the folder containing the FASTA files."
	exit 1
fi
input_folder="$1"

# Verify if directory is valid
if [ ! -d "$input_folder" ]; then
	echo "Folder '$input_folder' doesn't exist."
	exit 1
fi

# Name for the output
output="combined_sequences.fasta"

# Clean output if already exists
> $output

# Iterate through all the fastas found in path provided
for fasta_file in "$input_folder"/*.fasta; do
	# Obtain the name of the file without extension
	file_name=$(basename "$fasta_file" .fasta)
    
	# Process each fasta
	while IFS= read -r line; do
        # Modify the headers only
	if [[ $line == ">"* ]]; then
		# Add the name of the file as a first field in the header
		echo ">${file_name} ${line:1}" >> $output
	else
		# If is not a header, append without modifying
		echo "$line" >> $output
	fi
	done < "$fasta_file"
done
