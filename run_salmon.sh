#!/bin/bash

# Define directories
PAIRED_FASTQ_DIR="/path/paired_reads"  # Directory containing paired FASTQ files
SALMON_INDEX="/path/homosap_index"         # Path to Salmon index
SALMON_OUTPUT_BASE="/path/salmon_output"   # Base directory for Salmon output

# Create the base output directory if it doesn't exist
mkdir -p "$SALMON_OUTPUT_BASE"

# Loop through all _1_paired.fastq.gz files in the paired folder
for FILE1 in "$PAIRED_FASTQ_DIR"/*_1_paired.fastq.gz; do
    # Get the corresponding _2_paired file
    FILE2="${FILE1/_1_paired.fastq.gz/_2_paired.fastq.gz}"

    # Check if corresponding _2 file exists
    if [ -f "$FILE2" ]; then
        # Extract the base SRR ID (e.g., SRR14646315)
        BASE_NAME=$(basename "$FILE1" _1_paired.fastq.gz)

        # Define the output directory for this sample
        SAMPLE_OUTPUT_DIR="$SALMON_OUTPUT_BASE/$BASE_NAME"

        # Create the sample-specific output directory
        mkdir -p "$SAMPLE_OUTPUT_DIR"

        # Run Salmon quantification for this sample
        salmon quant -i "$SALMON_INDEX" \
            -l A \
            -1 "$FILE1" \
            -2 "$FILE2" \
            -p 8 \
            --validateMappings \
            -o "$SAMPLE_OUTPUT_DIR"

        echo "Salmon quantification completed for $BASE_NAME"
    else
        echo "Warning: No corresponding _2_paired.fastq.gz file found for $FILE1"
    fi
done