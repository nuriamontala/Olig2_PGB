#!/bin/bash

# This script downloads all the sequences (now it's fixed to proteins) provided in a list of accession numbers as first positional argument.
# It also names each file with the name list provided as second positional argument. Accession number and filename should be in the same position in the list, and both lists have the same length.

# Verify two arguments are passed.
if [ "$#" -ne 2 ]; then
	echo "Use: $0 <accessions_list> <names_list>"
	exit 1
fi

# Set variables
ACCESSION_FILE="$1"
NAMES_FILE="$2"


# Verify that both lists exist
if [ ! -f "$ACCESSION_FILE" ]; then
	echo "Error: Accession list '$ACCESSION_FILE' doesn't exist."
	exit 1
fi

if [ ! -f "$NAMES_FILE" ]; then
	echo "Error: Names list '$NAMES_FILE' doesn't exist."
	exit 1
fi

# Iterate through each row of the printed lists. While read allows to set the two columns used as input in different variables. 

paste "$ACCESSION_FILE" "$NAMES_FILE" | while IFS=$'\t' read -r accession name; do
# print message when processing each file
	echo "Downloading FASTA for Accession: $accession and saving as $name.fasta"
    
	# Download FASTA with efetch command and save with the desired filename. Set only to protein sequences.
	efetch -db protein -format fasta -id "$accession" > "${name}.fasta"
    
	# Verify if exit code is 0 meaning the sequence was succesfully downloaded
	if [ $? -eq 0 ]; then
		echo "File saved as ${name}.fasta"
	else
		echo "Unable to download ${accession}.fasta"
	fi
done

