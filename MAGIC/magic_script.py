import magic
import scprep

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

file = "counts_sorted.tab" # file in question

data = scprep.io.load_tsv(file,cell_axis='column',gene_names=True,cell_names=True)
# loading the tab file, cell_axis='column' as the format is genes x cells, gene_names & cell_names set to true as the file contains those
# data = scprep.filter.filter_empty_genes(data)
data = scprep.normalize.library_size_normalize(data) # library size normalization
data = scprep.transform.log(data, base = 2) # log transform with base 2

magic_op = magic.MAGIC() # we run MAGIC
magic_data = magic_op.fit_transform(data) # the new matrix of all genes (smoothed data

magic_pca = magic_op.transform(genes="pca_only") # extracting the principal components of the smoothed data
data_pca = PCA(n_components=3).fit_transform(np.array(data)) # PCA on raw data for comparison

# we plot principal components in 2d, before and after MAGIC
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,6))
scprep.plot.scatter2d(data_pca, label_prefix="PC", title='PCA without MAGIC', ax=ax1, ticks=False) # ta bort 'tick=False' om du vill se markeringar på koordinataxlarna
scprep.plot.scatter2d(magic_pca, label_prefix="PC", title='PCA with MAGIC', ax=ax2, ticks=False)
plt.tight_layout()
plt.show()

# we plot principal components in 3d, before and after MAGIC
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16, 6), subplot_kw={'projection':'3d'})
scprep.plot.scatter3d(data_pca, label_prefix="PC", title='PCA without MAGIC', ax=ax1, ticks=False)
scprep.plot.scatter3d(magic_pca, label_prefix="PC", title='PCA with MAGIC', ax=ax2, ticks=False)
plt.tight_layout()
plt.show()

magic_data = magic_op.transform(genes = "all_genes") # extract full smoothed matrix for downstream analysis
magic_data = magic_data.T # transpose data so that we have genes x cells
#magic_data.to_csv("counts_imputed.csv") # export data as a .csv file
