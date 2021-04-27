import magic
import scprep

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

file1 = 'counts_P0.tab' # file in question
file2 = 'counts_P7.tab'

data1 = scprep.io.load_tsv(file1,cell_axis='column',gene_names=True,cell_names=True) # use .load_csv for csv files
data2 = scprep.io.load_tsv(file2,cell_axis='column',gene_names=True,cell_names=True)
# loading the tab file, cell_axis='column' as the format is genes x cells, gene_names & cell_names set to true as the file contains those
# data = scprep.filter.filter_empty_genes(data)

A = scprep.stats.differential_expression(data1, data2, measure = 'ttest') # diff exp after importing count data

data1 = scprep.normalize.library_size_normalize(data1) # library size normalization
data2 = scprep.normalize.library_size_normalize(data2) # library size normalization
data1 = scprep.transform.log(data1, base = 2) # log transform with base 2
data2 = scprep.transform.log(data2, base = 2) # log transform with base 2

B = scprep.stats.differential_expression(data1, data2, measure = 'ttest') # diff exp after performing library size normalization and log transform (base 2)

magic_op_1 = magic.MAGIC() # we run MAGIC
magic_op_2 = magic.MAGIC()
magic_data1 = magic_op_1.fit_transform(data1) # the new matrix of all genes (smoothed data
magic_data2 = magic_op_2.fit_transform(data2) # the new matrix of all genes (smoothed data

C = scprep.stats.differential_expression(magic_data1, magic_data2, measure = 'ttest') # diff exp after performing MAGIC

print(A,B,C)
