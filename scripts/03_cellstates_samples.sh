#!/bin/bash

# Define the base directory and script path
BASE_DIR="/cellstatesappli/"
SCRIPT_PATH="/run_cellstates.py"

# List of samples
ITEMS=('Case1.YF' 'Case1.ZY' 'Case2.YF' 'Case2.ZC' 'Case2.ZY' 'Case3.YF' 'Case3.ZY' 'Case4.ZY')

# Loop through the samples
for item in "${ITEMS[@]}"; do

    # Define the input and output file paths
    INPUT_FILE="${BASE_DIR}RNAmatrix_${item}.tsv"
    OUTPUT_DIR="${BASE_DIR}RNAmatrix_${item}/"

    # Run the cellstates script
    python3 "$SCRIPT_PATH" "$INPUT_FILE" -t 18 -o "$OUTPUT_DIR"
done