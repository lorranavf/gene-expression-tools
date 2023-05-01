# Gene Expression Tools: Python functions for gene count  and statistics from RNA-Seq libraries downloaded from NCBI database
This repository contains a script in Python to perform RNA-Seq analysis on raw reads downloaded from the NCBI database. The script was developed using Python version 3.6.

### Function get_gene_count
Make sure you have installed all the dependencies and the required programs. The script expects the raw reads files to be downloaded from the NCBI database and saved in the readsRaw directory. The script also expects the adapter file to be saved in the adapters directory.
The script will perform the following steps:

- Download raw reads using fasterq-dump.
- Perform quality control analysis using FastQC.
- Perform adapter trimming using Trimmomatic.
- Perform quality control analysis on the trimmed reads using FastQC.
- Count reads per gene using Salmon.
- Generate general statistics using MultiQC.

The respective outputs of the analyses are in their respective directories.

The script also has a function named get_salmon_index that generates the index required for the Salmon program. This function takes as input the following 

### Requirements
Before running the script, make sure you have installed the following dependencies:

- Python version 3.6 or higher
- pandas version 1.2.0 or higher
- numpy version 1.19.5 or higher
- requests version 2.25.1 or higher
- seaborn version 0.11.1 or higher
- tqdm version 4.60.0 or higher
- matplotlib version 3.3.4 or higher

The following programs must also be installed and added to the system path:

- FastQC
- Trimmomatic
- fasterq-dump
- MultiQC
- Salmon
