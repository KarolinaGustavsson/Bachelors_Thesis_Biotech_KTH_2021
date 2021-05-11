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

components = 4                                              # choose nr of components in PCA
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

P0 = ['01', '02', '03', '04', '05', '06', '07', '08']
P7 = ['09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
red = ['01', '02', '03', '04', '09', '10', '11', '12', '13', '14', '15', '16']
green = ['05', '06', '07', '08', '17', '18', '19', '20', '21', '22', '23', '24']
P0_red = ['01', '02', '03', '04']
P0_green = ['05', '06', '07', '08']
P7_red = ['09', '10', '11', '12', '13', '14', '15', '16']
P7_green = ['17', '18', '19', '20', '21', '22', '23', '24']

# --------------------------------------------

def seperate_samples(samples, A, B):

    counts_list = samples.columns.tolist()                               # create list of column names in counts-file
    list_A = []
    list_B = []
    for x in counts_list:  # create lists for 2 column categories
        list_A.append(x if x[1:] in A else False)
        list_B.append(x if x[1:] in B else False)

    cols_A = samples[samples.columns.intersection(list_A)].columns   # select columns for 2 categories
    cols_B = samples[samples.columns.intersection(list_B)].columns

    samples_A = samples[cols_A]
    samples_B = samples[cols_B]

    return samples_A, samples_B

def create_new_df(table):
    cells = pd.DataFrame(table.columns.tolist())
    cells.columns = ['Cells']

    table = table.T
    table = table.reset_index(drop=True)
    table.columns = ["PC%s" %i for i in range(1,components+1)]
    table = cells.join(table)

    return table

def sep_component(sample, A, B):
    sample = sample.T                                           # transform sample matrix
    headers = sample.iloc[0]                                    # extract first row and store as headers
    sample  = pd.DataFrame(sample.values[1:], columns=headers)  # add back headers as headers for df

    sample_A, sample_B = seperate_samples(sample, A, B)         # seperate samples based on condition A and B
    sample_A = create_new_df(sample_A)                          # create df of sample A
    sample_B = create_new_df(sample_B)                          # create df of sample B

    # print(sample_A)
    # print('Amount of cells: ', len(sample_A.index))
    # print('               ')
    # print(sample_B)
    # print('Amount of cells: ', len(sample_B.index))
    # print('               ')

    return sample_A, sample_B

PCs = create_new_df(PCs)                                    # create df of components

# --------------------------------------------

# sort components and seperate into positive and negative component-value

sort_PC1 = PCs.sort_values(by=['PC1'], ascending=False)     # sort by PC1
sort_PC1 = sort_PC1.reset_index(drop=True)                  # reset index
positive_PC1 = sort_PC1.iloc[:97]                           # create df with positive values
negative_PC1 = sort_PC1.iloc[98:]                           # create df with negative values

sort_PC2 = PCs.sort_values(by=['PC2'], ascending=False)     # sort by PC2
sort_PC2 = sort_PC2.reset_index(drop=True)                  # reset index
positive_PC2 = sort_PC2.iloc[:211]                          # create df with positive values
negative_PC2 = sort_PC2.iloc[212:]                          # create df with negative values

sort_PC3 = PCs.sort_values(by=['PC3'], ascending=False)     # sort by PC3
sort_PC3 = sort_PC3.reset_index(drop=True)                  # reset index
positive_PC3 = sort_PC3.iloc[:274]                          # create df with positive values
negative_PC3 = sort_PC3.iloc[275:]                          # create df with negative values

sort_PC4 = PCs.sort_values(by=['PC4'], ascending=False)     # sort by PC4
sort_PC4 = sort_PC4.reset_index(drop=True)                  # reset index
positive_PC4 = sort_PC4.iloc[:142]                          # create df with positive values
negative_PC4 = sort_PC4.iloc[143:]                          # create df with negative values

# --------------------------------------------

# Check if component could be P0 vs P7

PC1_P0, PC1_P7 = sep_component(sort_PC1, P0, P7)
PC2_P0, PC2_P7 = sep_component(sort_PC2, P0, P7)
PC3_P0, PC3_P7 = sep_component(sort_PC3, P0, P7)
PC4_P0, PC4_P7 = sep_component(sort_PC4, P0, P7)

PC1_red, PC1_green = sep_component(sort_PC1, red, green)
PC2_red, PC2_green = sep_component(sort_PC2, red, green)
PC3_red, PC3_green = sep_component(sort_PC3, red, green)
PC4_red, PC4_green = sep_component(sort_PC4, red, green)

# --------------------------------------------

def plot_histogram(sample1, sample2, component, color1, color2, filename, bins):
    fig, ax = plt.subplots()
    ax.hist([sample1[component], sample2[component]], bins, color=[color1, color2], alpha=0.8)
    ax.set_ylabel('Number of cells')
    ax.set_xlabel('Component value')
    ax.grid(alpha=0.1)

    plt.suptitle('Number of cells in ' + component, fontsize='x-large', fontweight='medium')
    fig.tight_layout(pad=0.5, w_pad=2.0, h_pad=0.5)
    fig.legend(['P0 cells', 'P7 cells', 'P0 cells', 'P7 cells'], loc='upper right', bbox_to_anchor=(0.99, 0.92))
    plt.savefig(plots_dir + filename, dpi=200)
    plt.show()

def plot_histogram2(sample1, sample2, component, color1, color2, filename, bins):
    fig, ax = plt.subplots()
    ax.hist([sample1[component], sample2[component]], bins, color=[color1, color2], alpha=0.8)
    ax.set_ylabel('Number of cells')
    ax.set_xlabel('Component value')
    ax.grid(alpha=0.1)

    plt.suptitle('Number of cells in ' + component, fontsize='x-large', fontweight='medium')
    fig.tight_layout(pad=0.5, w_pad=2.0, h_pad=0.5)
    fig.legend(['Fucci red', 'Fucci green', 'Fucci red', 'Fucci green'], loc='upper right', bbox_to_anchor=(0.99, 0.92))
    plt.savefig(plots_dir + filename, dpi=200)
    plt.show()

plot_histogram(PC1_P0, PC1_P7, 'PC1', 'tab:pink', 'indigo', 'PC1_histogram.jpeg', np.linspace(-30, 60, 15))
plot_histogram(PC2_P0, PC2_P7, 'PC2', 'tab:pink', 'indigo', 'PC2_histogram.jpeg', np.linspace(-30, 30, 15))
plot_histogram(PC3_P0, PC3_P7, 'PC3', 'tab:pink', 'indigo', 'PC3_histogram.jpeg', np.linspace(-30, 30, 15))
plot_histogram(PC4_P0, PC4_P7, 'PC4', 'tab:pink', 'indigo', 'PC4_histogram.jpeg', np.linspace(-10, 10, 15))

plot_histogram2(PC1_red, PC1_green, 'PC1', 'tab:red', 'tab:green', 'PC1_fucci_histogram.jpeg', np.linspace(-30, 60, 15))
plot_histogram2(PC2_red, PC2_green, 'PC2', 'tab:red', 'tab:green', 'PC2_fucci_histogram.jpeg', np.linspace(-30, 30, 15))
plot_histogram2(PC3_red, PC3_green, 'PC3', 'tab:red', 'tab:green', 'PC3_fucci_histogram.jpeg', np.linspace(-30, 30, 15))
plot_histogram2(PC4_red, PC4_green, 'PC4', 'tab:red', 'tab:green', 'PC4_fucci_histogram.jpeg', np.linspace(-10, 10, 15))

# --------------------------------------------

# Remove high-valued PC1 cells and check PC2 for fucci again:

# high_PC1 = sort_PC1.iloc[:96]                   # seperate out high PC1
# high_PC1_cells = high_PC1['Cells']
# high_PC1_cells = np.asarray(high_PC1_cells)
#
# res = ''
# for ele in high_PC1_cells:
#     res = res + str(ele) + ","
#
# high_PC1_cells = res.split(",")
# high_PC1_cells.sort()
# high_PC1_cells.pop(0)
#
# counts = pd.read_csv('counts_imputed.csv', delimiter = ',', header=0, dtype='a')   # import data from csv-file
# genes = counts.iloc[:, :1]
# counts = counts.drop(['gene'], axis=1)
# cell_list = counts.columns.tolist()
#
# for element in high_PC1_cells:
#     if element in cell_list:
#         cell_list.remove(element)
#
# low_PC1_cells = cell_list
#
# high_PC1 = counts[counts.columns.intersection(high_PC1_cells)]
# low_PC1 = counts[counts.columns.intersection(low_PC1_cells)]
#
# high_PC1 = genes.join(high_PC1)
# low_PC1 = genes.join(low_PC1)
#
# high_PC1.to_csv('high_PC1.csv')
# low_PC1.to_csv('low_PC1')
#
# cond = PCs['Cells'].isin(high_PC1['Cells'])     # filter out cells from PCs
# PCs.drop(PCs[cond].index, inplace = True)       # remove high PC1 from PCs
# PCs = PCs.reset_index(drop=True)                # reset matrix index
#
# sort_PC2 = PCs.sort_values(by=['PC2'], ascending=False)     # sort by PC1
# sort_PC2 = sort_PC2.reset_index(drop=True)                  # reset index
#
# positive_PC2 = sort_PC2.iloc[:115]                          # create df with positive values
# negative_PC2 = sort_PC2.iloc[116:]                          # create df with negative values
#
# print('PC2:')
# sep_component(positive_PC2, red, green)
# sep_component(negative_PC2, red, green)

# --------------------------------------------

# Remove high-valued PC1 cells and check PC3 for P0 vs P7 again:

# high_PC1 = sort_PC1.iloc[:98]                   # seperate out high PC1
# cond = PCs['Cells'].isin(high_PC1['Cells'])     # filter out cells form PCs
# PCs.drop(PCs[cond].index, inplace = True)       # remove high PC1 from PCs
# PCs = PCs.reset_index(drop=True)                # reset matrix index
#
# sort_PC3 = PCs.sort_values(by=['PC3'], ascending=False)     # sort by PC3
# sort_PC3 = sort_PC3.reset_index(drop=True)                  # reset index
#
# positive_PC3 = sort_PC3.iloc[:177]                          # create df with positive values
# negative_PC3 = sort_PC3.iloc[178:]                          # create df with negative values
#
# print('PC3:')
# sep_component(positive_PC3, P0, P7)
# sep_component(negative_PC3, P0, P7)
