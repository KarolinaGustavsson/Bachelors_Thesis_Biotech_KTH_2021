import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os
script_dir = os.path.dirname("__file__")
plots_dir = os.path.join(script_dir, 'Plots/')
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)

# --------------------------------------------

file = 'imputed_data.h5ad'                                  # choose file to import
adata = scvi.data.read_h5ad(file)                           # import h5ad-file

components = 10                                             # choose nr of components in PCA
sc.pp.pca(adata, n_comps=components)                        # run PCA from scanpy

# --------------------------------------------

# Extract and prep principle components from AnnData

PCs = pd.DataFrame(adata.obsm['X_pca'])                     # extract PCA data from AnnData
PCs.columns = ["PC%s" %i for i in range(1,components+1)]    # rename columns according to principle components

cells = pd.DataFrame(adata.obs['cells'])                    # extract cell data from AnnData
cells.columns = ['Cells']                                   # rename DataFrame
cells = cells.reset_index(drop=True)                        # reset DataFrame index
PCs = PCs.reset_index(drop=True)                            # reset DataFrame index

PCs = cells.join(PCs).T                                     # join and transpose DataFrame
PCs = PCs.reset_index(drop=True)                            # reset DataFrame index

headers = PCs.iloc[0]                                       # extract top row
PCs  = pd.DataFrame(PCs.values[1:], columns=headers)        # turn top row into DataFrame header

# --------------------------------------------

# Create lists of cell groups

P0 = ['01', '02', '03', '04', '05', '06', '07', '08']
P7 = ['09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
red = ['01', '02', '03', '04', '09', '10', '11', '12', '13', '14', '15', '16']
green = ['05', '06', '07', '08', '17', '18', '19', '20', '21', '22', '23', '24']
P0_red = ['01', '02', '03', '04']
P0_green = ['05', '06', '07', '08']
P7_red = ['09', '10', '11', '12', '13', '14', '15', '16']
P7_green = ['17', '18', '19', '20', '21', '22', '23', '24']

# --------------------------------------------

# define functions

def seperate_samples(samples, A, B):

    counts_list = samples.columns.tolist()                   # create list of column names
    list_A = []
    list_B = []
    for x in counts_list:                                    # create lists for 2 column categories
        list_A.append(x if x[1:] in A else False)
        list_B.append(x if x[1:] in B else False)

    cols_A = samples[samples.columns.intersection(list_A)].columns   # select columns for 2 categories
    cols_B = samples[samples.columns.intersection(list_B)].columns

    samples_A = samples[cols_A]
    samples_B = samples[cols_B]

    return samples_A, samples_B

def create_new_df(table):
    cells = pd.DataFrame(table.columns.tolist())                # create list of column names
    cells.columns = ['Cells']                                   # add column title

    table = table.T                                             # transpose matrix
    table = table.reset_index(drop=True)                        # reset index
    table.columns = ["PC%s" %i for i in range(1,components+1)]  # name columns
    table = cells.join(table)                                   # join cell-name df with PC-data df

    return table

def sep_component(sample, A, B):
    sample = sample.T                                           # transform sample matrix
    headers = sample.iloc[0]                                    # extract first row and store as headers
    sample  = pd.DataFrame(sample.values[1:], columns=headers)  # add back headers as headers for df

    sample_A, sample_B = seperate_samples(sample, A, B)         # seperate samples based on condition A and B
    sample_A = create_new_df(sample_A)                          # create df of sample A
    sample_B = create_new_df(sample_B)                          # create df of sample B
    return sample_A, sample_B

def plot_PCA(A, B, filename):
    s = 12          # set marker size
    alpha = 0.75    # set marker opacity

    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(x = PCs_P0_red[A], y = PCs_P0_red[B], c='red', s=s, alpha=alpha)
    ax.scatter(x = PCs_P0_green[A], y = PCs_P0_green[B], c='green', s=s, alpha=alpha)
    ax.scatter(x = PCs_P7_red[A], y = PCs_P7_red[B], c='red', marker='v', s=s, alpha=alpha)
    ax.scatter(x = PCs_P7_green[A], y = PCs_P7_green[B], c='green', marker='v', s=s, alpha=alpha)
    ax.set_xlabel(A)
    ax.set_ylabel(B)
    ax.grid(alpha=0.7)
    ax.legend(['P0 (red)', 'P0 (green)', 'P7 (red)', 'P7 (green)'], loc='lower right', bbox_to_anchor=(0.99, 0.01))

    fig.tight_layout(pad=1.5, w_pad=2, h_pad=1)
    plt.savefig(plots_dir + filename, dpi=200)
    plt.show()

# --------------------------------------------

# Remove positive PC1 cells

# PCs = create_new_df(PCs)
# sort_PC1 = PCs.sort_values(by=['PC1'], ascending=False)     # sort by PC1
# sort_PC1 = sort_PC1.reset_index(drop=True)                  # reset index
# high_PC1 = sort_PC1.iloc[:98]                               # seperate out high PC1
# cond = PCs['Cells'].isin(high_PC1['Cells'])                 # filter out cells form PCs
# PCs.drop(PCs[cond].index, inplace = True)                   # remove high PC1 from PCs
# PCs = PCs.reset_index(drop=True)                            # reset DataFrame index
# PCs = PCs.T                                                 # transpose DataFrame
# PCs = PCs.reset_index(drop=True)                            # reset DataFrame index
# headers = PCs.iloc[0]                                       # extract top row
# PCs  = pd.DataFrame(PCs.values[1:], columns=headers)        # turn top row into DataFrame header

# --------------------------------------------

# Remove positive PC4 cells

# PCs = create_new_df(PCs)
# sort_PC4 = PCs.sort_values(by=['PC4'], ascending=False)     # sort by PC1
# sort_PC4 = sort_PC4.reset_index(drop=True)                  # reset index
# high_PC4 = sort_PC4.iloc[:142]                              # seperate out high PC1
# cond = PCs['Cells'].isin(high_PC4['Cells'])                 # filter out cells form PCs
# PCs.drop(PCs[cond].index, inplace = True)                   # remove high PC1 from PCs
# PCs = PCs.reset_index(drop=True)                            # reset matrix index
# PCs = PCs.T                                                 # transpose DataFrame
# PCs = PCs.reset_index(drop=True)                            # reset DataFrame index
# headers = PCs.iloc[0]                                       # extract top row
# PCs  = pd.DataFrame(PCs.values[1:], columns=headers)        # turn top row into DataFrame header

# --------------------------------------------

# seperate and create new DataFrames for the samples

PCs_P0, PCs_P7 = seperate_samples(PCs, P0, P7)
PCs_P0_red, PCs_P0_green = seperate_samples(PCs_P0, P0_red, P0_green)
PCs_P7_red, PCs_P7_green = seperate_samples(PCs_P7, P7_red, P7_green)

PCs_P0_red = create_new_df(PCs_P0_red)
PCs_P0_green = create_new_df(PCs_P0_green)
PCs_P7_red = create_new_df(PCs_P7_red)
PCs_P7_green = create_new_df(PCs_P7_green)
PCs = create_new_df(PCs)

# --------------------------------------------

# Plot the first four components

plot_PCA('PC1', 'PC2', 'PC1_PC2.jpeg')
plot_PCA('PC1', 'PC3', 'PC1_PC3.jpeg')
plot_PCA('PC2', 'PC3', 'PC2_PC3.jpeg')
plot_PCA('PC4', 'PC2', 'PC4_PC2.jpeg')

# --------------------------------------------
