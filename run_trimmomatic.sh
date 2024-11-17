#!/bin/bash

# Define directories
FASTQ_DIR="/path/fastqs"
OUTPUT_DIR="/path/outputs"
TRIMMOMATIC_JAR="/path/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTER_FILE="/path/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all _1.fastq.gz files in the fastqs directory
for FILE1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    # Get the corresponding _2 file
    FILE2="${FILE1/_1.fastq.gz/_2.fastq.gz}"

    # Check if corresponding _2 file exists
    if [ -f "$FILE2" ]; then
        # Get the base name (e.g., SRR14646315)
        BASE_NAME=$(basename "$FILE1" _1.fastq.gz)

        # Define output file names
        OUTPUT_1_PAIRED="$OUTPUT_DIR/${BASE_NAME}_1_paired.fastq.gz"
        OUTPUT_1_UNPAIRED="$OUTPUT_DIR/${BASE_NAME}_1_unpaired.fastq.gz"
        OUTPUT_2_PAIRED="$OUTPUT_DIR/${BASE_NAME}_2_paired.fastq.gz"
        OUTPUT_2_UNPAIRED="$OUTPUT_DIR/${BASE_NAME}_2_unpaired.fastq.gz"
        LOG_FILE="$OUTPUT_DIR/${BASE_NAME}_trimlog.txt"

        # Run Trimmomatic for the pair
        java -jar "$TRIMMOMATIC_JAR" PE -threads 4 -phred33 \
            -trimlog "$LOG_FILE" \
            "$FILE1" "$FILE2" \
            "$OUTPUT_1_PAIRED" "$OUTPUT_1_UNPAIRED" \
            "$OUTPUT_2_PAIRED" "$OUTPUT_2_UNPAIRED" \
            ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
            SLIDINGWINDOW:5:20 \
            MINLEN:50
    else
        echo "Warning: No corresponding _2.fastq.gz file found for $FILE1"
    fi
done