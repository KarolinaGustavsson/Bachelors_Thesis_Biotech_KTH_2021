import pandas as pd
import numpy as np
import sklearn

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




# more attempts at splitting
# lines = table
# for line in lines:
#     line.slice(',')[16].split('\t')

# print( list(x for x in table.strip().split('\t'))  for x in table)

# table2 = table.strip("\t")
# table3 = table.reader(table, delimiter="\t")
# print(table3.head())

# for aRow in table:
#     print aRow.split('\t')
# we tried
# tabstuff = ["a"	"b"	"c"	"d"	"e"	"f"	"g",
# 4	65	87	98	56	43	2,
# 45	67	67	8	90	34	23,
# 23	43	43	56	78	89	98]
# print(tabstuff)