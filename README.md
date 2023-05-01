# Gene Expression Tools: Python functions for quantify the expression of transcripts  and statistics from RNA-Seq libraries downloaded from NCBI database
This repository contains a script in Python to perform RNA-Seq analysis on raw reads downloaded from the NCBI database. The script was developed using Python version 3.6.

### Function: get_salmon_index

This function is a Python function used to create an index for Salmon. The function takes in the paths to the genome and transcriptome fasta files, a gtf file, and an optional gff file. It first converts the gff file to gtf using the gffread command if the gff file is provided. Then it generates a decoy transcriptome using the **generateDecoyTranscriptome.sh** script from Salmon repository. Finally, it creates an index for Salmon using the salmon index command with the decoy transcriptome fasta file and the decoys.txt file generated by generateDecoyTranscriptome.sh.

To call this function:

```python
from rna_seq_pipeline import get_salmon_index
genome = '/path/to/genome.fa'
transcriptome = '/path/to/transcriptome.fa'
generate_decoy_script= '/path/to/generateDecoyTranscriptome.sh'
salmon = '/path/to/Salmon'
gtf_file = '/path/to/gtf_file.gtf'
get_salmon_index(genome, transcriptome, generate_decoy_script, salmon,gtf_file)
```

### Function: get_gene_count
Make sure you have installed all the dependencies and the required programs. The script expects the raw reads files to be downloaded from the NCBI database and saved in the readsRaw directory. The script also expects the adapter file to be saved in the adapters directory and the index directory.
The script will perform the following steps:

- Download raw reads using fasterq-dump.
- Perform quality control analysis using FastQC.
- Perform adapter trimming using Trimmomatic.
- Perform quality control analysis on the trimmed reads using FastQC.
- Count reads per gene using Salmon.
- Generate general statistics using MultiQC.

The respective outputs of the analyses are in their respective directories.

To call this function:
```python
from rna_seq_pipeline import get_gene_count
workdir = '/path/to/workdir'
samples = ['sample1', 'sample2', 'sample3']
adapters = '/path/to/adapters.fa'
index = '/path/to/index'
fasterq = '/path/to/fasterq-dump'
fastqc = '/path/to/FastQC'
trimmomatic = '/path/to/Trimmomatic.jar'
salmon = '/path/to/Salmon'
multiqc = '/path/to/MultiQC'
get_gene_count(workdir, samples, adapters, index, fasterq, fastqc, trimmomatic, salmon, multiqc)
```

### Function: get_table_count

This function takes a directory path as input, containing Salmon output files, and returns a pandas DataFrame with TPM values for each sample. It first gets a list of samples and folders containing the output files, and then checks if the quant.sf file is present in each folder. It then reads in each table and filters the columns to include only the 'Name' and 'TPM' columns, and then renames the 'TPM' column with the corresponding sample name. Finally, it concatenates all the tables and returns the TPM values for each sample as a DataFrame.

To call this function:
```python
from rna_seq_pipeline import get_table_count
df_count = get_table_count("path/to/salmon/output/directory")
```

### Function: get_corr_var_to_var

This function takes a table with numerical values, a correlation method ('pearson', 'kendall', or 'spearman'), and a correlation coefficient threshold as input. It then calculates the correlation between variables (columns) in the table and returns a DataFrame with the variables that have a correlation coefficient above the given threshold. First, it removes columns that have a sum of zero, then calculates the correlation matrix for the remaining columns using the specified method. It then filters the matrix to include only the upper triangle (excluding the diagonal) and stores the names of the columns with a correlation coefficient above the threshold. Finally, it returns a DataFrame with the names of the two variables and their correlation coefficient.

To call this function:
```python
from rna_seq_pipeline import get_corr_var_to_var
corr_var = get_corr_var_to_var(df_count, "pearson", 0.8)
```

### Function: generate_heatmap

This function takes a pandas DataFrame as input and generates a heatmap using the Seaborn library. It first creates a figure and axis object with the specified size, then generates the heatmap using the 'crest' colormap. It then sets the x and y axis labels and the title of the plot, and finally displays the plot using plt.show().

To call this function:
```python
from rna_seq_pipeline import generate_heatmap
generate_heatmap(df_count, xlabel="Samples", ylabel="Genes", title="TPM Heatmap")
```

### Function: generate_clustermap

This function takes a pandas DataFrame as input and generates a clustermap using the Seaborn library. It first generates a clustermap using the 'crest' colormap and specifying the position and size of the colorbar. It then displays the plot using plt.show().

To call this function:
```python
from rna_seq_pipeline import generate_clustermap
generate_clustermap(df_count)
```

### Requirements
Before running the script, make sure you have installed the following dependencies:

- Python version 3.6 or higher
- pandas version 1.2.0 or higher
- numpy version 1.19.5 or higher
- requests version 2.25.1 or higher
- tqdm version 4.60.0 or higher
- matplotlib version 3.3.4 or higher
- seaborn version 0.11.1 or higher

The following programs must also be installed and added to the system path:

- [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools/wiki)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [MultiQC](https://multiqc.info/)
- [Salmon](https://combine-lab.github.io/salmon/)
- [GffRead](https://github.com/gpertea/gffread)

generateDecoyTrnascriptome.sh dependencies:
- awk
- [bedtools](https://github.com/arq5x/bedtools2)
- mashmap
