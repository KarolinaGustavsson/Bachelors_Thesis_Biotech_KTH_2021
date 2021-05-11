import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------

pd.set_option('display.max_rows', None)                 # show unlimited amount of rows in output window

# --------------------------------------------

file = 'imputed_data.h5ad'                              # choose file to import
adata = scvi.data.read_h5ad(file)                       # import h5ad-file

components = 10                                         # choose nr of components in PCA
sc.pp.pca(adata, n_comps=components)                    # run PCA from scanpy

# --------------------------------------------

# Extract variance from anndata

variance_ratio = pd.DataFrame(adata.uns['pca']['variance_ratio']).T                         # extract variance ratio from AnnData
variance_ratio.columns = ["PC%s" %i for i in range(1,components+1)]                         # rename columns according to principle component
variance_ratio_PC = pd.DataFrame(variance_ratio.columns.tolist())                           # save column names
variance_ratio_PC.columns = ['Components']                                                  # rename DataFrame
variance_ratio = variance_ratio.T                                                           # transpose DataFrame
variance_ratio = variance_ratio.reset_index(drop=True)                                      # reset DataFrame index
variance_ratio.columns = ['Variance ratio']                                                 # rename DataFrame
variance_ratio = variance_ratio_PC.join(variance_ratio)                                     # join variance ratios with component names
cumulative_variance_ratio = pd.DataFrame(np.cumsum(variance_ratio['Variance ratio']))       # calculate cumulative variance
cumulative_variance_ratio.columns = ['Cumulative variance']                                 # rename DataFrame
variance = variance_ratio.join(cumulative_variance_ratio)                                   # join all components

# --------------------------------------------

# Plot the figure as bar and line chart

fig, ax = plt.subplots()
ax.bar(variance['Components'], variance['Variance ratio'])
ax.plot(variance['Components'], variance['Cumulative variance'], 'r-')
ax.set_ylabel('Explained variance')
ax.set_ylim([0.0, 1.1])
ax.grid(alpha=0.1)

plt.suptitle('Explained varience per principle component', fontsize='x-large', fontweight='medium')
fig.tight_layout(pad=0.5, w_pad=2.0, h_pad=0.5)
fig.legend(['Cumulative variance', 'Explained variance'], loc='upper right', bbox_to_anchor=(0.99, 0.8))
plt.savefig('PCA_explained_variance.jpeg', dpi=200)
plt.show()
