import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy

table = pd.read_csv("counts.tab")
# table = pd.read_csv("magictabs.rtf")
print(table.head())
print(np.shape(table))
print(np.ndim(table))
print(np.prod(table.shape))
length = np.prod(table.shape)
print(length)

data1 = pd.DataFrame(table)
print(data1.head)

# line = list(range(0, int(length), 1))
# print(line)

# this one below works! but only for the first row...
# not useful, only parses first row, too complicated
for row in table:
    print( list(zip(*(line.strip().split('\t') for line in table))) )

print("first try")
# omg omg omg this works!
# omg omg omg this works!
# delimiter = tabs, separates each value
table4 = pd.read_csv("counts.tab", delimiter = '\t')
print(table4.head())
print("parsed")
# [5 rows x 385 columns]
# actual size is [24490 rows x 385 columns]
table5 = table4.T
table6 = table5.drop('gene') #tar bort översta raden
table7 = table6.values
print(table7)
print("Yay it works!")

#THIS ACTUALLY WORKS, YEET - BUT WE SHOULD TRY ANOTHER METHOD AS WELL
normalized_X = preprocessing.normalize(table7)
print(normalized_X)
standardized_X = preprocessing.scale(table7)
print(standardized_X)

print("LMAO PLZ WORK")
#OTHER METHOD FOR STANDARDIZATION
import scipy
from scipy.sparse import csr_matrix
table8 = scipy.sparse.csr_matrix(table7)
arr = np.array(table8)
print(csr_matrix(arr))

#print("plz work")
#normalized_X = preprocessing.normalize(table10)
#print(normalized_X)
# funkade inte
#table7 = scipy.sparse.csr_matrix(table6.values)
#print(table7.head())

# import scipy
# from scipy.sparse import csr_matrix
# arr = np.array(table5)
# print(csr_matrix(arr))

# we want to convert it into a sparse matrix
# scaled_data = preprocessing.scale(table5)
