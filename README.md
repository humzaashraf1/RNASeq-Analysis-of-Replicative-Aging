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

<img src="https://github.com/user-attachments/assets/afa60f81-15d7-43e7-8812-a4a0c4182b61" alt="fastqc_per_base_sequence" width="400"/>

<img src="https://github.com/user-attachments/assets/f25deace-2cef-4ec6-89b1-3d752fe7fab4" alt="fastqc_per_base_sequence" width="430"/>

The overrepresented sequences tab suggests that the TruSeq universal Illumina adapters are causing significant read contamination. To deal with this we can use **Trimmomatic** (http://www.usadellab.org/cms/?page=trimmomatic), a java-based adapter trimming tool that comes pre-packaged with Illumina TruSeq.fa files. Since I am on windows, I am using Windows Subsystem for Linux (https://learn.microsoft.com/en-us/windows/wsl/install) to run this command.
