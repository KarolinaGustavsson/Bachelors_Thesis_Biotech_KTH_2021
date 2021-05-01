import time                                         #import timer
startTime = time.time()                             #timer start

# --------------------------------------------

import pandas as pd                                 # import pandas
import numpy as np                                  # import numpy
import matplotlib.pyplot as plt                     # import matplotlib (for plotting)
from sklearn import preprocessing                   # import part of sklearn
from sklearn.decomposition import PCA               # import PCA from sklearn
import seaborn as sns                               # import seaborn
from numpy import inf                               # import part of numpy
import numpy.random as npr                          # import random from numpy
from scipy.stats import ttest_ind                   # import ttest for hypothesis testing

# --------------------------------------------

# This part imports the data from CSV-file

file = "counts_fucci.tab"                                           # choose what file to import data from

counts = pd.read_csv(file, delimiter = '\t', header=0, dtype='a')   # import data from csv-file
counts = counts.drop(['gene'], axis=1)                              # remove first column (gene)

# --------------------------------------------

# To remove fucci-info use this part!

table_fucci = counts.iloc[:1, :]                    # create matrix with only fucci-data
table_fucci = table_fucci.T                         # transpose the matrix
table_fucci = table_fucci.reset_index(drop=True)    # reset dataframe index

counts = counts.iloc[1: , :]                          # remove fucci-data from matrix
counts = counts.apply(pd.to_numeric)

# counts.loc["A0"] = counts.loc['01']

counts = counts.reset_index(drop=True)                # reset dataframe index
counts.dropna(axis=0, how='all', inplace=True)
counts = counts.loc[~(counts<=0.0).all(axis=1)]
counts = pd.DataFrame(data=np.log2(counts),index=counts.index,columns=counts.columns)
counts[counts == -inf] = 0

print(counts)


#
# # --------------------------------------------

# This part creates three dataframes, counts (with all data), counts_P0 (with P0 data) and counts_P7 (with P7 data)
counts_list = counts.columns.tolist()                               # create list of column names in counts-file
list_cols_1_8 = []
list_cols_9_x = []
for x in counts_list:  # create lists for 2 column categories
    list_cols_1_8.append(x if x[1:] in ['01', '02', '03', '04', '05', '06', '07', '08'] else False)
    list_cols_9_x.append(x if x[1:] not in ['01', '02', '03', '04', '05', '06', '07', '08'] else False)

cols_1_8 = counts[counts.columns.intersection(list_cols_1_8)].columns   # select columns for 2 categories
cols_9_x = counts[counts.columns.intersection(list_cols_9_x)].columns


def get_significance_two_groups(row):
    log_fold_change = row[cols_1_8].mean() - row[cols_9_x].mean()  # Calculate the log Fold Change
    p = ttest_ind(row[cols_1_8], row[cols_9_x], equal_var=False)[1]  # Calculate the significance # False --> Welsh
    return [p, -np.log10(p), log_fold_change]


pvalues = counts.apply(get_significance_two_groups, axis=1, result_type="expand")
pvalues.rename(columns={list(pvalues)[0]: 'p', list(pvalues)[1]: '-log_p', list(pvalues)[2]: 'log_FC'}, inplace=True)

print(pvalues)

# plt.figure(figsize=(12, 8))
# sns.displot(pvalues["p"],kde=False) # plot distribution of p-values
# plt.xlim(0,1.0);
# plt.xlabel("P-Value")
# plt.ylabel("Cell Count")
# plt.title("P-Value Distribution")


# sns.set_style("white")
# sns.set_context("talk")
# ax = sns.relplot(data=pvalues,x="log_FC",y="-log_p",aspect=1.5,height=6)
# # plt.axhline(p_treshold)
# #sns.lineplot([-6,4],[p_treshold,p_treshold],ax=ax)
# ax.set(xlabel="$log_2(TN/not TN)$", ylabel="$-log_{10}(p)$");
# plt.show()

def bootstrap(invec):
    idx = npr.randint(0, len(invec), len(invec))
    return [invec[i] for i in idx]

def estimatePi0(p, numBoot=100, numLambda=100, maxLambda=0.95):
    p.sort()
    n=len(p)
    lambdas=np.linspace(maxLambda/numLambda,maxLambda,numLambda)
    Wls=np.array([n-np.argmax(p>=l) for l in lambdas])
    pi0s=np.array([Wls[i] / (n * (1 - lambdas[i])) for i in range(numLambda)])
    minPi0=np.min(pi0s)
    mse = np.zeros(numLambda)
    for boot in range(numBoot):
        pBoot = bootstrap(p)
        pBoot.sort()
        WlsBoot =np.array([n-np.argmax(pBoot>=l) for l in lambdas])
        pi0sBoot =np.array([WlsBoot[i] / (n *(1 - lambdas[i])) for i in range(numLambda)])
        mse = mse + np.square(pi0sBoot-minPi0)
    minIx = np.argmin(mse)
    return pi0s[minIx]

def qvalues(pvalues):
    m = pvalues.shape[0]  # The number of p-values
    pvalues.sort_values("p", inplace=True)  # sort the p-values in ascending order
    pi0 = estimatePi0(list(pvalues["p"].values))
    print("pi_0 estimated to " + str(pi0))

    # calculate a FDR(t) as in Storey & Tibshirani
    num_p = 0.0
    for ix in pvalues.index:
        num_p += 1.0
        t = pvalues.loc[ix, "p"]
        fdr = pi0 * t * m / num_p
        pvalues.loc[ix, "q"] = fdr
        pi0_hat = (m - num_p) / (m * (1 - t))
        pvalues.loc[ix, "pi0_hat"] = pi0_hat

    # calculate a q(p) as the minimal FDR(t)
    old_q = 1.0
    for ix in reversed(list(pvalues.index)):
        q = min(old_q, pvalues.loc[ix, "q"])
        old_q = q
        pvalues.loc[ix, "q"] = q
    return pvalues

qv = qvalues(pvalues) #vi får ut ungefär 0.88 -> 12% av generna är differentially expressed

print("vi vill toppgenerna:")
print(qv.iloc[0:40,0:4]) #40 är antal rader vi vill se

sns.lineplot(pvalues["q"],list(range(pvalues.shape[0])),ci=None,lw=3)
plt.xlim(0,0.1);
plt.ylim();
plt.ylabel("Number of differential genes");
plt.show()

#qv #det här tar jättelång tid, oklart om det funkar, ska ge signifikans för alla gener

qv["Significant"] = qv["q"]<1e-10
less_than_FDR_10 = qv[qv["q"]<1e-10]
p_treshold = float(less_than_FDR_10.iloc[-1:]["-log_p"].values)

sns.set_style("white")
sns.set_context("talk")
ax = sns.relplot(data=pvalues,x="log_FC",y="-log_p",hue="Significant",aspect=1.5,height=6)
plt.axhline(p_treshold)
#sns.lineplot([-6,4],[p_treshold,p_treshold],ax=ax)
ax.set(xlabel="$log_2(TN/not TN)$", ylabel="$-log_{10}(p)$");
plt.show()

# # --------------------------------------------

executionTime = (time.time() - startTime)/60                                       #timer stop.
print('Execution time in seconds: ' + str(executionTime))

