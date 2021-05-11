
# --------------------------------------------

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# --------------------------------------------

file = 'counts_imputed.csv'                                        # choose file to import
data = pd.read_csv(file, delimiter = ',', header=None, dtype='a')  # import data as dataframe with float64

# --------------------------------------------

data.drop([0], axis=1, inplace=True)                            # remove first column
obs = data.iloc[:1, :]                                          # extract cells as variable obs
obs = obs.drop([1], 1).T                                        # remove first column and transpose obs
var = data.iloc[:, :1]                                          # extract genes as variable var
var = var.drop(0, axis=0)                                       # remove first column of var
var.reset_index(drop=True, inplace=True)                        # reset obs index
obs.reset_index(drop=True, inplace=True)                        # reset var index

# --------------------------------------------

data.drop([1], axis=1, inplace=True)                            # remove first column
data.drop([0], axis=0, inplace=True)                            # remove first row
data.reset_index(drop=True, inplace=True)                       # reset data index
data = data.T                                                   # transpose dataframe
data.reset_index(drop=True, inplace=True)                       # reset data index

# --------------------------------------------

obs.columns = ['cells']                                         # rename observation "cells"
var.columns = ['genes']                                         # rename variables "cells"
data.to_numpy()                                                 # convert df "data" to numpy array

# --------------------------------------------

adata = ad.AnnData(data, obs=obs, var=var, dtype='float64')     # create AnnData object
adata.write('top200.h5ad')                                      # create h5ad-file

# --------------------------------------------
