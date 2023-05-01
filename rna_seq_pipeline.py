#!/usr/bin/env python

import os
import time
import pandas as pd
import numpy as np
import requests
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List, Union

def get_salmon_index(genome:str, transcriptome:str, generate_decoy_script:str, salmon:str, gtf_file:str, gff_file:str=None) -> None:
    """
    Generates an index for Salmon.

    Args:
        genome (str): Path to the genome file.
        transcriptome (str): Path to the transcriptome file.
        gtf_file (str): Path to the gtf file.
        gff_file (str): Path to the gff file. Default is None.
        generate_decoy_script (str): Path to the script for generating decoy sequences.
        salmon (str): Path to the Salmon executable.

    Returns:
        None.
    
    Obs:
        Generates an index for Salmon on index directory.
    """

    # get gtf
    if gff_file:
        os.system(f'gffread {gff_file} -T -o {gtf_file}')
    else:
        pass

    # get decoy
    os.system(f'bash {generate_decoy_script} -a {gtf_file} -g {genome} -t {transcriptome}  -o decoy')

    # get index
    os.system(f'{salmon} index -t decoy/gentrome.fa -i index --decoys decoy/decoys.txt -k 31')


def get_gene_count(workdir: str, samples:List[str], adapters:str, index:str, fasterq:str, fastqc:str, trimmomatic:str, salmon:str, multiqc:str) -> None:
    """
    Quantify trnascripts and general statistics from RNA-Seq libraries downloaded from NCBI database using the specified programs.

    Args:
        workdir (str): Path to the working directory.
        samples (list): List of sample names.
        adapters (str): Path to the adapter Trimmomatic file.
        index (str): Path to the index Salmon directory.
        fasterq (str): Path to the fasterq-dump program.
        fastqc (str): Path to the FastQC program.
        trimmomatic (str): Path to the Trimmomatic program.
        salmon (str): Path to the Salmon program.
        multiqc (str): Path to the MultiQC program.

    Returns:
        None. 
        
    Obs:
        The respective outputs of the analyses are in their respective directories.
    """
    
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    else:
        pass

    # Reads Raw
    raw_left = [workdir + '/readsRaw/' + i + '_1.fastq' for i in samples]
    raw_right = [workdir + '/readsRaw/' + i + '_2.fastq' for i in samples]

    # Reads Clean
    clean_left = [workdir + '/readsClean/' + i + '_1_pair.fq' for i in samples]
    clean_right = [workdir + '/readsClean/' + i + '_2_pair.fq' for i in samples]
    clean_left_unpair = [workdir + '/readsClean/' + i + '_1_unpair.fq' for i in samples]
    clean_right_unpair = [workdir + '/readsClean/' + i + '_2_unpair.fq' for i in samples]

    # Diretórios das análises
    os.system(f'mkdir {workdir}/FastQC_raw')
    os.system(f'mkdir {workdir}/FastQC_clean')
    os.system(f'mkdir {workdir}/MultiQC_raw')
    os.system(f'mkdir {workdir}/MultiQC_clean')
    os.system(f'mkdir {workdir}/MultiQC_salmon')
    os.system(f'mkdir {workdir}/Salmon')

    def check_internet():
        url = 'https://www.google.com'
        timeout = 5
        try:
            requests.get(url, timeout=timeout)
            return True
        except requests.exceptions.ConnectionError:
            return False

    for i in range(len(samples)):

        # check se internet on
        internet = 'offline'

        while internet == 'offline':
            if not check_internet():
                internet = 'offline'
                time.sleep(60)
            else:
                internet = 'online'
            if internet == 'online':
                break

        # make dir
        os.system(f'mkdir {workdir}/readsRaw')
        os.system(f'mkdir {workdir}/readsClean')
        os.system(f'mkdir {workdir}/Salmon/{samples[i]}')

        # start log
        log1 = 'Iniciando análise da amostra: {} [{}]'.format(i + 1, samples[i])

        # fasterq-dump log
        log2 = 'Fasterq-dump iniciou às     : {}'.format(time.ctime())

        # fasterq-dump analysis
        os.system(f'{fasterq} {samples[i]} -m 10000MB -e 8 -t {workdir}/readsRaw --outdir {workdir}/readsRaw')

        # fastqc log
        log3 = 'FastQC raw iniciou às       : {}'.format(time.ctime())

        # fastqc analysis
        os.system(f'{fastqc} {raw_left[i]}  -t 8 -o {workdir}/FastQC_raw')
        os.system(f'{fastqc} {raw_right[i]} -t 8 -o {workdir}/FastQC_raw')

        # trimmomatic log
        log4 = 'Trimmomatic iniciou às      : {}'.format(time.ctime())

        # trimmomatic analysis
        os.system(f'java -jar {trimmomatic} PE -threads 4 -phred33 {raw_left[i]} {raw_right[i]} {clean_left[i]} {clean_left_unpair[i]} {clean_right[i]} {clean_right_unpair[i]} ILLUMINACLIP:{adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')

        # fastqc log
        log5 = 'FastQC clean inicou às      : {}'.format(time.ctime())

        # fastqc analysis
        os.system(f'{fastqc} {clean_left[i]}  -t 8 -o {workdir}/FastQC_clean')
        os.system(f'{fastqc} {clean_right[i]} -t 8 -o {workdir}/FastQC_clean')

        # salmon log
        log6 = 'Salmon iniciou às           : {}'.format(time.ctime())

        # salmon analysis
        os.system(f'{salmon} quant -i {index} -l IU -1 {clean_left[i]} -2 {clean_right[i]} -o {workdir}/Salmon/{samples[i]} -q')

        # remove sample dir
        os.system(f'rm -r {workdir}/readsRaw/')
        os.system(f'rm -r {workdir}/readsClean/')

        # end log
        log7 = 'Encerrando as análises às   : {}'.format(time.ctime())

        # save logs
        logs = open(workdir + '/log.out', 'a+')
        logs.write('\n'.join([log1, log2, log3, log4, log5, log6, log7]) + '\n')
        logs.close()

    # multiqc log
    log9 = 'MultiQC iniciou às          : {}'.format(time.ctime())

    # multiqc analysis
    os.system(f'{multiqc} {workdir}/FastQC_raw/. --outdir {workdir}/MultiQC_raw')
    os.system(f'{multiqc} {workdir}/FastQC_clean/. --outdir {workdir}/MultiQC_clean')
    os.system(f'{multiqc} {workdir}/Salmon/. --outdir {workdir}/MultiQC_salmon')

    log10 = 'MultiQC encerrou às         : {}'.format(time.ctime())

    # end analysis
    log11 = 'Análises finalizadas às     : {}'.format(time.ctime())

    # save logs
    logs = open(workdir + '/log.out', 'a+')
    logs.write('\n'.join([log9, log10, log11]) + '\n')
    logs.close()


def get_table_count(salmon_outdir:str) -> pd.DataFrame:
    """
    Given a list of samples in a directory containing salmon output files, returns a DataFrame with TPM values
    for each sample.

    Args:
        salmon_outdir (str): Directory path containing salmon output files.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the TPM values for each sample.
    """


    ext = '/quant.sf'
    samples = os.listdir(salmon_outdir)
    folders = [f'{salmon_outdir}/{i}' for i in samples]
    tables = [f'{salmon_outdir}/{i}/{ext}' for i in samples]

    error = []

    for i in range(len(samples)):

        ls = os.listdir(path=folders[i])

        if 'quant.sf' not in ls:
            error.append(samples[i])
            print('Repeat analysis for sample {}'.format(samples[i]))
        else:
            pass

    dfs = []

    for i in range(len(samples)):

        if samples[i] not in error:
            df = pd.read_table(tables[i])
            df = df.filter(items=['Name', 'TPM'])
            df = df.rename(columns={'TPM': samples[i]})
            df = df.set_index(['Name'])
            dfs.append(df)

        else:
            pass

    df_count = pd.concat(dfs, axis=1, join='inner').T

    return df_count


def get_corr_var_to_var(table:Union[str, pd.DataFrame], method:str, treshold:float) -> pd.DataFrame:
    """
    Calculates the correlation between variables (columns) in a table and returns a dataframe with the 
    variables (columns) that have a correlation coefficient above a given threshold.

    Args:
    - table: str or pd.DataFrame - A table with numerical values.
    - method: str - A string indicating the correlation method. Options are 'pearson', 'kendall' and 'spearman'.
    - treshold: float - A float between 0 and 1 indicating the minimum correlation coefficient to be included.

    Returns:
    - A pandas dataframe with the following columns:
        - Var1: str - Name of the first variable in the correlation pair.
        - Var2: str - Name of the second variable in the correlation pair.
        - Correlation: float - Correlation coefficient between Var1 and Var2.
    """


    if isinstance(table, pd.DataFrame):
        df = table.copy()
    else:
        df = pd.read_csv(table)

    # remove gene if gene sum is zero
    df = df.iloc[:, list(df.sum() != 0)]

    # get corr table for genes
    dfc = df.corr(method=method)
    dfc = dfc.to_numpy(na_value=0)
    dfc = np.triu(dfc, 1)

    # get gene gene interaction best correlation
    var1 = []
    var2 = []
    corr = []
    genes = list(df.columns)

    for i in tqdm(range(len(dfc))):
        for j in range(len(dfc)):
            if abs(dfc[i, j]) >= treshold:
                var1.append(genes[i])
                var2.append(genes[j])
                corr.append(dfc[i, j])

    corr_var_to_var = pd.DataFrame({'Var1': var1, 'Var2': var2, 'Correlation': corr})

    return corr_var_to_var

def generate_heatmap(df:pd.DataFrame, xlabel:str=None, ylabel:str=None, title:str=None) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(df, cmap='crest')
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    plt.show()
    
def generate_clustermap(df:pd.DataFrame) -> None:
    sns.clustermap(df, cmap="crest", cbar_pos=(1, 0.15, .025, .5))
    plt.show()
