# conda create -n dgi_coexpre python=3.9
import pandas as pd
import numpy as np
# import numba
# from numba import jit
import sys, os
sys.path.append('..')
import pickle as pkl

# import data
from data import load_counts_groups, Load_OutputData_df, DESeq, Get_DESeq_label, obtain_idname
from config import config

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# load read counts and groups data for preprocessing, finally obtain count dataframe and group information (conditions)
counts_df, conditions = load_counts_groups(config)

# DESeq analysis 
dds, stat_res = DESeq(config, counts_df, conditions, SAVE=1) 

# load RPKM data (normalized gene abundance)
RPKM_matrix = pd.read_table(config['RPKM_PATH'], sep='\t', header=0, index_col=0)

# labeling DEGs
df = Get_DESeq_label(conditions, counts_df, config, RPKM_matrix, # MATRIX_COL_NAME, 
                    dds, SAVE=1)
