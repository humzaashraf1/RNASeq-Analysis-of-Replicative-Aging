# RNASeq Analysis of Replicative Aging
In this tutorial, we will explore differential gene expression in immortalized versus non-immortalized WI-38 fibroblasts at varying population doublings. The replicative lifespan of mammalian cells is limited by the end-replication problem, a gradual loss of genomic material from chromosome ends during cell division. Chromosome ends are capped by telomeres—repeated "TTAGGG" sequences that shorten with age. When telomeres become critically short, cells enter an irreversible, non-dividing state called senescence, a process associated with many aging-related diseases. This tutorial aims to distinguish molecular features that arise due to proliferation alone from those caused by replicative aging (telomere shortening). To do so, we’ll compare h-TERT (telomerase reverse transcriptase) "immortalized" primary cells, which do not lose telomere length with population doublings, to replicatively aged cells after the same number of cell divsions. This will allow us to differentiate molecular effects due to cell division alone versus replicative decay, deepening our understanding of cellular aging.

The dataset used here comes from an excellent paper published in eLife by researchers at Calico called **Novel insights from a multiomics dissection of the Hayflick limit** (https://doi.org/10.7554/eLife.70283). We will go through step-by-step how to process and analyze some of the data shown in Figure 1E. The raw data for the paper is available through SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP321317. The GEO Accession link with the processed data is here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175533.

To start, we will use **SRAToolKit** (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) to download 12 samples using the **fastq-dump** command. The easiest way to do this is to create fetch_seqs.bat and execute it using the command line. The script requires you to specify the location of fastq-dump.exe, the path to an output folder for downloading the fastq.gz files, and a .txt file containing only the SRA IDs with no header for each sample:

| Sample Name       | SRA ID       |
|-------------------|--------------|
| hTERT_TP1_A       | SRR14646263  |
| hTERT_TP1_B       | SRR14646264  |
| hTERT_TP1_C       | SRR14646265  |
| hTERT_TP5_A       | SRR14646272  |
| hTERT_TP5_B       | SRR14646273  |
| hTERT_TP5_C       | SRR14646274  |
| RS_PDL20_TP1_A    | SRR14646293  |
| RS_PDL20_TP1_B    | SRR14646294  |
| RS_PDL20_TP1_C    | SRR14646295  |
| RS_PDL50_TP8_A    | SRR14646311  |
| RS_PDL50_TP8_B    | SRR14646312  |
| RS_PDL50_TP8_C    | SRR14646313  |

fetch_seqs.bat:
```batch
@echo off
set "TOOL_PATH=C:\path\SRAToolKit\bin\fastq-dump.exe"
set "OUT_DIR=C:\path\Downloaded"
set "SRR_LIST=srr_list.txt"

if not exist "%TOOL_PATH%" (
    echo Error: fastq-dump not found at %TOOL_PATH%
    exit /b
)

if not exist "%SRR_LIST%" (
    echo Error: srr_list.txt not found in the current directory
    exit /b
)

for /F "tokens=*" %%s in (%SRR_LIST%) do (
    echo Downloading %%s...
    "%TOOL_PATH%" --split-files --gzip --outdir "%OUT_DIR%" %%s
    if %errorlevel% neq 0 (
        echo Failed to download %%s
    )
)
echo All downloads complete.
```
Next, we will run **fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to determine which adapters are present from the Illumina sequencing. After downloading the release, I unzipped the root and executed the run_fastqc.bat file. This prompted me to load in an unzipped fastq file, creating a quality control report in the GUI:

<img src="https://github.com/user-attachments/assets/afa60f81-15d7-43e7-8812-a4a0c4182b61" alt="fastqc_per_base_sequence" height = "300" width="400"/>

<img src="https://github.com/user-attachments/assets/f25deace-2cef-4ec6-89b1-3d752fe7fab4" alt="fastqc_per_base_sequence" height = "300" width="400"/>

The overrepresented sequences tab suggests that the TruSeq universal Illumina adapters are causing significant read contamination. To deal with this we can use **Trimmomatic** (http://www.usadellab.org/cms/?page=trimmomatic), a java-based adapter trimming tool that comes pre-packaged with Illumina TruSeq.fa files. Since I am on windows, I am using Windows Subsystem for Linux (https://learn.microsoft.com/en-us/windows/wsl/install) to run this command.

**From the Trimmomatic User Guide:**

<img src="https://github.com/user-attachments/assets/fc0d8000-2952-4ee0-9fef-6edecb10ab8a" alt="fastqc_per_base_sequence" height = "300" width="400"/>

"For paired-end data, two input files, and 4 output files are specified, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not."

"java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> | <unpaired output 1> <paired output 2> <unpaired output 2> <step 1> ..."

**Trimmomatic Parameters:**

ILLUMINACLIP:"$ADAPTER_FILE":2:30:10
Removes adapter sequences based on the provided adapter file. Allows up to 2 mismatches, with trimming occurring if the adapter has a quality score of 30 or above.

SLIDINGWINDOW:5:20
Trims bases from the read if the average quality score in a sliding window of 5 bases is below 20.

MINLEN:50
Discards reads that are shorter than 50 bases after trimming.

run_trimmomatic.sh:
```bash
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


```
Forward reads before versus after running Trimmomatic:

<img src="https://github.com/user-attachments/assets/f41e2a9b-0888-4c00-bfb5-24ecb786a12e" alt="fastqc_per_base_sequence" height = "300" width="400"/>
<img src="https://github.com/user-attachments/assets/61b91373-32d2-4f8a-b956-8ac9dca02589" alt="fastqc_per_base_sequence" height = "300"  width="400"/>

Reverse reads before versus after running Trimmomatic:

<img src="https://github.com/user-attachments/assets/be3b858e-a99b-4c3c-86af-6c236d4ee0f7" alt="fastqc_per_base_sequence" height = "300" width="400"/>
<img src="https://github.com/user-attachments/assets/1162a7de-cf2e-43d9-93ca-76fb8685423d" alt="fastqc_per_base_sequence" height = "300" width="400"/>

Adapter content before versus after running Trimmomatic:

<img src="https://github.com/user-attachments/assets/ee26dda5-41a5-4f98-b443-b71d9dd1bacb" alt="fastqc_per_base_sequence" height = "285" width="400"/>
<img src="https://github.com/user-attachments/assets/38d28309-c9d2-470b-b40d-053efdf677f3" alt="fastqc_per_base_sequence" height = "285" width="400"/>

<sub> Portions of code in this repository were generated with the assistance of ChatGPT, a LLM developed by OpenAI.</sub>
