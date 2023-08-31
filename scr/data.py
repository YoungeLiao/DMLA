#####pretreatment#####
# packages needed
import pandas as pd
import numpy as np
import numba
from numba import jit
import sys, os
sys.path.append('..')
import pickle as pkl

# bioinformatics
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import config

## hyperparameter/variables
SAVE = 1 # 1 represents True, indicating save the results; 0 represents False, indicating don't save results
# load data & group information

# functions

def load_counts_groups(config, filter_mean_thre=0):
    r"""
    load external data, please set the parameters at config.py
    --- input ---
    :para config['DATAPATH']: pathway of raw reads counts datasets. NOTE! We highly recommend organize the data matrix by put the reference (control) group on the first column
    :para config['FEATUREPATH']: pathway of features (conditions) datasets
    :para config['TEMP_PATH']: pathway of temperary datasets
    --- output ---
    :para config['counts_df']: file name of counts matrix
    :para conditions: condition and group information, i.e. features
    """
    # Read the table from a file
    data_path = config['DATAPATH']
    count_matrix = pd.read_table(data_path, sep='\t', header=0, index_col=0)
    counts_df = count_matrix.transpose()
    
    # set the filtering threshold
    if filter_mean_thre == 0:
        thre = 1/len(counts_df.iloc[:,0])
    else: 
        thre = filter_mean_thre

    # load group/condition information
    conditions = pd.read_csv(config['FEATUREPATH'], index_col=0)

    # filter out genes that have less than 10 read counts 
    genes_to_keep = counts_df.columns[counts_df.mean(axis=0) >= thre]
    counts_df = counts_df[genes_to_keep]
    
    # make output directory 
    if not os.path.exists(config['OUTPUT_PATH']):
        os.makedirs(config['OUTPUT_PATH'])
    config['TEMP_PATH'] =  config['OUTPUT_PATH'] + config['TEMP_FILENAME']
    
    # save data
    counts_df.to_csv(config['TEMP_PATH'], index=True)

    # return valid data & group information
    return counts_df, conditions

def obtain_idname(config):
    """
    :return idname: name of id, which could be utilized for merge dataframe by column name
    :return id: id values, i.e. the names of sequences
    """
    count_matrix = pd.read_table(config['DATAPATH'], sep='\t', header=0)
    idname = count_matrix.iloc[:,0].name
    id = count_matrix.iloc[:,0].values
    return idname, id

# Some useful documentation
# [annData](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
# [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/api/docstrings/pydeseq2.ds.DeseqStats.html)
def DESeq(config, counts_df, conditions, SAVE):
    r"""
    DESeq analysis to obtain basic results
    --- input ---
    :para counts_df: counts matrix
    :para conditions: condition and group information matrix, i.e. features
    :para config['factors']: belong to the conditions, e.g. group
    :para SAVE: whether or not save the returned variables as files, default is true
    :para config['OUTPUT_PATH']: directory for output data 
    --- output ---
    :para dds: DESeq results, raw results that need further statistical analysis
    :para stat_res: statistical analyisis on DESeq results across all groups
    """
    factor = config['factors']
    dds = DeseqDataSet(
        counts=counts_df,
        clinical=conditions,
        design_factors=factor,
        refit_cooks=True
    )
    
    # run the deseq2() method to fit dispersions and LFCs.
    dds.deseq2()
    
    # obtain p-value & fc: Statistical analysis with the DeseqStats class
    stat_res = DeseqStats(dds)
    stat_res.summary()
    if not os.path.exists(config['OUTPUT_PATH']):
        os.makedirs(config['OUTPUT_PATH'])
    if SAVE:
        with open(os.path.join(config['OUTPUT_PATH'], "dds.pkl"), "wb") as f:
            pkl.dump(dds, f)
        with open(os.path.join(config['OUTPUT_PATH'], "stat_results.pkl"), "wb") as f:
            pkl.dump(stat_res, f)

    return dds, stat_res

def DESeq_res_stat_paried(dds, counts_df, config, Test, SAVE):
    r"""
    DESeq analysis to obtain basic results
    --- input ---
    :para counts_df: counts matrix
    :para config['condition']: condition that utilized for statistical analysis, such as p-value calculation. For example, we utilized 'light' () to indicate grouped by illunination condition for statistical analysis 
    :para config['factors']: belong to the conditions, e.g. group
    :para SAVE: whether or not save the returned variables as files, default is true
    :para config['OUTPUT_PATH']: directory for output data 
    --- output ---
    :para stat_res: statistical results of paired statistical analysis
    :para pvalue_df: p-value results of the statistical results
    :para pvalue_filtered: only genes with p-value <= 0.05 were kept
    """
    stat_res = DeseqStats(dds, contrast=[config['condition'], Test, config['ref']])
    stat_res.summary()
    pvalue_df = pd.DataFrame(stat_res.p_values)
    pvalue_filtered = pvalue_df[pvalue_df.sum(axis=1) <= 0.05]
    if not os.path.exists(config['OUTPUT_PATH']):
        os.makedirs(config['OUTPUT_PATH'])
    if SAVE:
        with open(config['OUTPUT_PATH'] + Test +'_vs_' + config['ref'] + "_pvalue_filtered.pkl", 'wb') as f:
            pkl.dump(pvalue_filtered, f)
        with open(config['OUTPUT_PATH'] + Test +'_vs_' + config['ref'] + "_stat_results.pkl", "wb") as f:
            pkl.dump(stat_res, f)
    return stat_res, pvalue_df, pvalue_filtered


def Get_DESeq_label(conditions, counts_df, config, RPKM_matrix, 
                    dds, SAVE, filter_mean_thre=0):
    r"""
    DESeq analysis to obtain basic results
    --- input ---
    :para conditions: condition and group information matrix, i.e. features
    :para config['ref']: config['ref']erence group for DESeq analysis to obtain p-value
    :para RPKM_PATH: the normalized gene abundance matrix quantified by RPKM
    :para MATRIX_COL_NAME: column names of RPKM matrix
    :para dds: DESeq results
    :para config['condition']: condition that utilized for statistical analysis, such as p-value calculation. For example, we utilized 'light' () to indicate grouped by illunination condition for statistical analysis 
    --- output ---
    :para df: data matrix
    """
    groups = conditions.iloc[:, 0].unique()
    test_groups = np.delete(groups, np.where(groups == config['ref']))

    # set the filtering threshold
    if filter_mean_thre == 0:
        thre = 1/len(counts_df.iloc[:,0])
    else: 
        thre = filter_mean_thre

    # RPKM
    RPKM_matrix_filtered = RPKM_matrix[RPKM_matrix.mean(axis=1) >= thre]
    MATRIX_COL_NAME = RPKM_matrix_filtered.columns.tolist()
    idname, id = obtain_idname(config)
    df = pd.DataFrame(columns= [idname,'p_value', 'labels'] + MATRIX_COL_NAME) # initialize dataframe for cotaining all labeled data

    for test in test_groups:
        DESeq_res_stat_paried(dds, counts_df, config, test, SAVE)
        PVALUE_PKL_PATH = config['OUTPUT_PATH'] + test +'_vs_' + config['ref'] + "_pvalue_filtered.pkl"
        CONDITION_COL_NAME = [idname,'p_value']
        genes_p_filtered = Load_OutputData_df(PVALUE_PKL_PATH, CONDITION_COL_NAME)
        genes_p_filtered['labels'] = test + 'vs' + config['ref']
        df1 = genes_p_filtered
        df2 = RPKM_matrix_filtered
        intersected_df = pd.merge(df1, df2, on = idname)
        # intersected_df = pd.merge(df1, df2, left_index=True, right_index=True, how='inner')
        # df = df.append(intersected_df, ignore_index=True)
        df = pd.concat([df, intersected_df])
     
    if not os.path.exists(config['OUTPUT_PATH'] + 'ToRdata/'):
        os.makedirs(config['OUTPUT_PATH'] + 'ToRdata/')
    if SAVE:
        df.to_csv(config['OUTPUT_PATH'] + 'ToRdata/labeled_matrix.csv', index=False)

    return df

def Load_OutputData_df(file_path, COL_NAME, TYPE='pkl'):
    r"""
    DESeq analysis to obtain basic results
    --- input ---
    :para file_path: pathway of file to be loaded
    --- output ---
    :para data: data in dataframe
    """
    if TYPE == 'pkl':
        with open(file_path, 'rb') as file:
        # Load the contents of the .pkl file
            data = pkl.load(file)
            data.reset_index(inplace=True)
            data.columns = COL_NAME
    elif TYPE == 'csv':
        # Load the contents of the .csv file
        data = pd.read_csv(file_path)
    elif TYPE == 'txt':
        data = pd.read_table(file_path, sep='\t', header=0)
    return data


def get_existing_dataset(config, data_type='csv'):
    file_path = config['data_path']
    if data_type == 'pkl':
        with open(file_path, 'rb') as file:
        # Load the contents of the .pkl file
            data = pkl.load(file)
            data.reset_index(inplace=True)
    elif data_type == 'csv':
        # Load the contents of the .csv file
        data = pd.read_csv(file_path)
    elif data_type == 'txt':
        data = pd.read_table(file_path, sep='\t', header=0)
    return data


def FilterDesq(df, group, pvalue, fc, x_thre):
    r"""
    DESeq analysis to obtain basic results
    --- input ---
    :para df: labeled deseq results, i.e. the DEGs (differentially expressed genes) of different groups
    :para group: targeted filtering groups
    :para pvalue: filtering threshold of p value 
    :para fc: filtering threshold of fold changes
    :para dds: DESeq results
    :para config['condition']: condition that utilized for statistical analysis, such as p-value calculation. For example, we utilized 'light' () to indicate grouped by illunination condition for statistical analysis 
    --- output ---
    :para df: data matrix
    """
    a = 1
    return a

