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
After downloading, it is a good practice to clear your cache to not run into future issues with SRA-queries. To do this, navigate to your /path/SRAToolKit/bin folder and run the following command:

```
vdb-config
```

Next, we will run **fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to determine which adapters are present from the Illumina sequencing. After downloading the release, I unzipped the root and executed the run_fastqc.bat file. This prompted me to load in an unzipped fastq file, creating a quality control report in the GUI:

<img src="https://github.com/user-attachments/assets/afa60f81-15d7-43e7-8812-a4a0c4182b61" alt="1" height = "300" width="400"/>

<img src="https://github.com/user-attachments/assets/f25deace-2cef-4ec6-89b1-3d752fe7fab4" alt="2" height = "300" width="400"/>

The overrepresented sequences tab suggests that the TruSeq universal Illumina adapters are causing significant read contamination. To deal with this we can use **Trimmomatic** (http://www.usadellab.org/cms/?page=trimmomatic), a java-based adapter trimming tool that comes pre-packaged with Illumina TruSeq.fa files. Since I am on windows, I am using Windows Subsystem for Linux (https://learn.microsoft.com/en-us/windows/wsl/install) to run this command.

**From the Trimmomatic User Guide:**

<img src="https://github.com/user-attachments/assets/fc0d8000-2952-4ee0-9fef-6edecb10ab8a" alt="3" height = "300" width="400"/>

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

<img src="https://github.com/user-attachments/assets/f41e2a9b-0888-4c00-bfb5-24ecb786a12e" alt="4" height = "300" width="400"/>
<img src="https://github.com/user-attachments/assets/61b91373-32d2-4f8a-b956-8ac9dca02589" alt="5" height = "300"  width="400"/>

Reverse reads before versus after running Trimmomatic:

<img src="https://github.com/user-attachments/assets/be3b858e-a99b-4c3c-86af-6c236d4ee0f7" alt="6" height = "300" width="400"/>
<img src="https://github.com/user-attachments/assets/1162a7de-cf2e-43d9-93ca-76fb8685423d" alt="7" height = "300" width="400"/>

Adapter content before versus after running Trimmomatic:

<img src="https://github.com/user-attachments/assets/ee26dda5-41a5-4f98-b443-b71d9dd1bacb" alt="8" height = "285" width="400"/>
<img src="https://github.com/user-attachments/assets/38d28309-c9d2-470b-b40d-053efdf677f3" alt="9" height = "285" width="400"/>

After trimming our reads, we are ready to map our reads and generate a counts matrix. An efficient way to do this for gene expression is through an open-source tool called Salmon (https://combine-lab.github.io/salmon/). Salmon performs quasi-mapping directly on the transcriptome, greatly reducing the computational resources required to process our data. The best way to access Salmon is through Bioconda on a Linux OS. There is a good tutorial in their documentation for getting the enviornment set up. The first step is to download the human transcriptome and create an index (in a directory of your choice):

```bash
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz -o homosap.fa.gz
```
```bash
salmon index -t homosap.fa.gz -i homosap_index
```
After creating our Conda-env, we can call on Salmon by passing our index file, forward and reverse reads, and a directory to write the output:
```bash
salmon quant -i /path/homosap_index \
-l A \
-1 /path/outputs/SRRID_1_paired.fastq.gz \
-2 /path/outputs/SRRID_2_paired.fastq.gz \
-p 8 \
--validateMappings \
-o /path/salmon_output
```
Here is a bash script to achieve the same result for a folder of paired reads (run_salmon.sh):
```bash
#!/bin/bash

# Define directories
PAIRED_FASTQ_DIR="/path/to/paired_folder"  # Directory containing paired FASTQ files
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
```
Next, we can compile our quant.sf outputs for each sample using **pytximport** (https://pytximport.readthedocs.io/en/dev/start.html). I created a conda enviornment with WSL and wrote a small script in python to aggregate and export the raw data to a single csv file for all the samples (get_expression.py):
```python
from pytximport import tximport
from pytximport.utils import create_transcript_gene_map
import pandas as pd
import os

quant_files = []
salmon_path = '/path/salmon_output'
for folder in os.listdir(salmon_path):
    quant_path = salmon_path + '/' + folder + '/' + 'quant.sf'
    quant_files.append(quant_path)

transcript_gene_map = create_transcript_gene_map(species="human")

results = tximport(
    file_paths=quant_files,
    data_type="salmon",
    transcript_gene_map=transcript_gene_map
)

expression_df = pd.DataFrame(
    results.X.toarray() if hasattr(results.X, "toarray") else results.X,
    index=results.obs.index, 
    columns=results.var.index
)
combined_df = pd.concat([results.obs, expression_df], axis=1)
combined_df.to_csv("counts_matrix.csv")
```
The output should look something like this:
|                 | ENSG00000123457 | ENSG00000123458 | ENSG00000123459 | ENSG00000123460 | ... |
|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|
| Sample1            | 12.56           | 99.67           | 5.23            | 23.45           | ...          |
| Sample2            | 8.43            | 87.34           | 3.67            | 19.23           | ...           |
| Sample3            | 1.23            | 56.89           | 8.67            | 34.56           | ...           |
| Sample4            | 3.45            | 67.89           | 2.34            | 0.00            | ...            |
| ...                | 0.00            | 45.67           | 6.78            | 22.33           | ...            |

We can use this table as an input for differential gene-expression analysis to **pydeseq2** (https://pydeseq2.readthedocs.io/en/latest/), a python implementation of DeSeq2 in R. Differential gene expression analysis involves comparing sample-specific mean expression levels to global mean expression levels across multiple replicates. This approach identifies genes that are statistically significantly enriched or depleted under different biological conditions. The paper describing the method can be found here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
<img src="https://github.com/user-attachments/assets/09d28e56-f887-410e-8c42-2257c0531e9f" alt="10" height = "300" width="600"/>

run_deseq2.ipynb:
```python
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

counts_df = pd.read_csv('/path/counts_matrix.csv')
counts_df.rename(columns={'Unnamed: 0': 'SampleID'}, inplace=True)

data = {
    "Sample Name": [
        "hTERT_TP1", "hTERT_TP1", "hTERT_TP1",
        "hTERT_TP5", "hTERT_TP5", "hTERT_TP5",
        "RS_PDL20_TP1", "RS_PDL20_TP1", "RS_PDL20_TP1",
        "RS_PDL50_TP8", "RS_PDL50_TP8", "RS_PDL50_TP8"
    ],
    "SRA ID": [
        "SRR14646263", "SRR14646264", "SRR14646265",
        "SRR14646272", "SRR14646273", "SRR14646274",
        "SRR14646293", "SRR14646294", "SRR14646295",
        "SRR14646311", "SRR14646312", "SRR14646313"
    ]
}
data_df = pd.DataFrame(data)

counts_df["SRA ID"] = counts_df["SampleID"].str.extract(r'/([^/]+)/quant\.sf$')[0]

metadata = counts_df.merge(data_df, on="SRA ID", how="left")
metadata = metadata[["SRA ID", "Sample Name"]].rename(columns={"Sample Name": "condition"})

gene_names = pd.DataFrame(counts_df.columns, columns=["Gene Names"])
gene_names = gene_names[~gene_names["Gene Names"].isin(["SRA ID", "SampleID"])]

counts_df = counts_df.drop(columns=["SRA ID", "SampleID"])
counts_df.columns = range(counts_df.shape[1])
counts_df = counts_df.round().astype(int)
counts_df = counts_df.apply(pd.to_numeric, errors="raise")
```
With the data correctly formatted and the relevant metadata extracted, we can set up the deseq model:
```python
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
    inference=inference,
)

dds.deseq2()

dds.obs
```
To get the statistically significant differential gene expression for a pair of samples:
```python
stat_res = DeseqStats(dds, contrast = ('condition','hTERT-TP1','hTERT-TP5'))
stat_res.summary()
res = stat_res.results_df
res['ensembl'] = gene_names['Gene Names'].values

sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
sigs
```
<sub> Portions of code in this repository were generated with the assistance of ChatGPT, a LLM developed by OpenAI.</sub>
