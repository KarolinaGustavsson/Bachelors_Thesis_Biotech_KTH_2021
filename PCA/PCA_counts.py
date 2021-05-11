# -----------------------------------------------

import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

import plotly.graph_objs as go
import plotly.offline as offline
offline.init_notebook_mode()

# -----------------------------------------------

table_counts = pd.read_csv("counts.tab", delimiter = '\t')
#print(table_counts.head())

#print("-------------------------")

table_counts = table_counts.T
table_counts = table_counts.drop('gene') #tar bort översta raden
table_counts = table_counts.values
table_counts = table_counts.T
#print(table_counts)

print("-------------------------")

normalized_counts = preprocessing.normalize(table_counts)

#print("-------------------------")

components = 11
pca = PCA(n_components=components)
Y = pca.fit_transform(normalized_counts)
#var_exp = Y.explained_variance_ratio_  #ignore
#cum_var_exp = np.cumsum(var_exp)       #ignore

# Save components to a DataFrame
PCA_components = pd.DataFrame(Y)

#print(PCA_components)

# --------------------------------------------

# Plot the explained variances

# features = range(pca.n_components_)
# plt.bar(features, pca.explained_variance_ratio_, color='black')
# plt.xlabel('PCA features')
# plt.ylabel('variance %')
# plt.xticks(features)

# --------------------------------------------

# Använd nedanstående för att plot i en 2D graf

plt.scatter(PCA_components[0], PCA_components[1], alpha=.1, color='black')
plt.xlabel('PCA 1')
plt.ylabel('PCA 2')

# --------------------------------------------


# Allt nedanför är extra

# # Plot the explained variance
# x = ["PC%s" %i for i in range(1,components)]
# trace1 = go.Bar(
#     x=x,
#     y=list(var_exp),
#     name="Explained Variance")
#
# trace2 = go.Scatter(
#     x=x,
#     y=cum_var_exp,
#     name="Cumulative Variance")
#
# layout = go.Layout(
#     title='Explained variance',
#     xaxis=dict(title='Principle Components', tickmode='linear'))
#
# data = [trace1, trace2]
# fig = go.Figure(data=data, layout=layout)
# fig
#
# # ----------------------------------------
#
# # Project first three components
# # Y_train_pca = pca.fit_transform(normalized_counts)
# #
# # traces = []
# # for name in ['ALL', 'AML']:
# #     trace = go.Scatter3d(
# #         x=Y_train_pca[y_train.cancer==name,0],
# #         y=Y_train_pca[y_train.cancer==name,1],
# #         z=Y_train_pca[y_train.cancer==name,2],
# #         mode='markers',
# #         name=name,
# #         marker=go.Marker(size=10, line=go.Line(width=1),opacity=1))
# #
# #     traces.append(trace)
# #
# # layout = go.Layout(
# #     xaxis=dict(title='PC1'),
# #     yaxis=dict(title='PC2'),
# #     title="Projection of First Three Principle Components"
# # )
# #
# # data = traces
# # fig = go.Figure(data=data, layout=layout)
# # fig
