#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import plotly.figure_factory as ff
import plotly.graph_objs as go
import plotly.express as px
from plotly.offline import plot
from scipy.cluster import hierarchy as hc
from scipy.cluster.hierarchy import linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
import argparse as ap 

# from django.conf import settings


parser = ap.ArgumentParser(description="preps RNAseq counts table for analysis")
parser.add_argument('a', help="input file of RNAseq Counts")
args=parser.parse_args()


def processCounts(df):
    df = df.loc[~(df == 0).all(axis=1)]
    df2 = df.filter(regex='^(?!MIR)', axis=0)
    df3 = df2.filter(regex='^(?!__)', axis=0)
    return df3


# Function to Quantile normalise the counts matrix
def quantileNormalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)



# Log transform the the data
def log2Transform(df):
    df = df + 1
    df_log = df.applymap(np.log2)
    return df_log


pca_3d = open("pca_3d_rnaseq.html", 'w')
pca_1_2_2d = open("pca_1_2_2d_rnaseq.html", 'w')
pca_1_3_2d = open("pca_1_3_2d_rnaseq.html", 'w')
pca_1_4_2d = open("pca_1_4_2d_rnaseq.html", 'w')

def scatterPlot3D():

    # read in txt file of counts
    counts = pd.read_csv(args.a, index_col="Gene", sep="\t")

    counts = processCounts(counts)

    counts = quantileNormalize(counts)

    counts = log2Transform(counts)

    counts_t = counts.T

    counts_t_new = counts_t.index.str.split(" ", n=1, expand=True)

    

    kayrotype = []
    for i in counts_t_new:
        kayrotype.append(i[1])

    gene_name = []
    for i in counts_t_new:
        gene_name.append(i[0])



    counts_t["Gene"] = gene_name
    counts_t["Kayrotype"] = kayrotype


    X = counts_t.iloc[:, 0:-2].values
    y = counts_t.iloc[:, -1].values


    colors = {'MLL': '#F6EE03',
              'Complex': '#F67403',
              '-7/(del(7)': '#C60006',
              'T(6:9)': '#0006C6',
              'Others': '#FCADFF',
              '-5/del(5)': '#5D0060'}

    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=3)


    Y_sklearn = sklearn_pca.fit_transform(X_std)

    data = []

    for name, col in zip(('MLL', 'Complex', '-7/(del(7)', 'T(6:9)', 'Others', '-5/del(5)'), colors.values()):
        trace = dict(
            type='scatter3d',
            x=Y_sklearn[y == name, 0],
            y=Y_sklearn[y == name, 1],
            z=Y_sklearn[y == name, 2],
            mode='markers',
            name=name,
          #  hovertext=,
            marker=dict(
                color=col,
                size=12,
                line=dict(
                    color='rgba(0, 0, 0, 0.0)',
                    width=0.7),
                opacity=0.7)
        )
        data.append(trace)


    layout = {
        "plot_bgcolor": "rgba(255,35,35,1)",
        "scene": {
            "xaxis": {
                "type": "linear",
                "title": "PC1",

                
            },
            "yaxis": {
                "type": "linear",
                "title": "PC2"
            },
            "zaxis": {
                "type": "linear",
                "title": "PC3",
            },
        }
    }

    fig = go.Figure(data=data, layout=layout)

    fig.update_layout({'width': 800, 'height': 800,
                       'showlegend': True, 'hovermode': 'closest'},
                       margin=dict(l=50,r=50,b=100,t=100,pad=4))

    #fig.update_yaxes(automargin=True)

    pca3d_div = plot(fig, include_plotlyjs=True, output_type='div')

    return pca3d_div



def scatterPlot2D():

    # read in txt file of counts
    counts = pd.read_csv(args.a, index_col="Gene", sep="\t")

    counts = processCounts(counts)

    counts = quantileNormalize(counts)

    counts = log2Transform(counts)

    counts_t = counts.T

    counts_t_new = counts_t.index.str.split(" ", n=1, expand=True)

    kayrotype = []
    for i in counts_t_new:
        kayrotype.append(i[1])

    gene_name = []
    for i in counts_t_new:
        gene_name.append(i[0])

    counts_t["Gene"] = gene_name
    counts_t["Kayrotype"] = kayrotype

    X = counts_t.iloc[:, 0:-2].values
    y = counts_t.iloc[:, -1].values


    colors = {'MLL': '#F6EE03',
              'Complex': '#F67403',
              '-7/(del(7)': '#C60006',
              'T(6:9)': '#0006C6',
              'Others': '#FCADFF',
              '-5/del(5)': '#5D0060'}

    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=3)

    Y_sklearn = sklearn_pca.fit_transform(X_std)

    data = []

    for name, col in zip(('MLL', 'Complex', '-7/(del(7)', 'T(6:9)', 'Others', '-5/del(5)'), colors.values()):


        trace = dict(
                type='scatter',
                x=Y_sklearn[y == name, 0],
                y=Y_sklearn[y == name, 1],
                    # z=Y_sklearn[y == name, 2],
                mode='markers',
                name=name,
                marker=dict(
                    color=col,
                    size=12,
                    line=dict(
                        color='rgba(0, 0, 0, 0.0)',
                        width=0.7),
                    opacity=0.7)
                )
        data.append(trace)

        layout = dict(
                xaxis=dict(title='PC1', showline=False),
                yaxis=dict(title='PC2', showline=False)
        )


    fig = go.Figure(data=data, layout=layout)

    fig.update_layout({'width': 600, 'height': 600,
                       'showlegend': True, 'hovermode': 'closest'})

    pca2d_div = plot(fig, include_plotlyjs=True, output_type='div')

    return pca2d_div



def scatterPlot2D_1_3():

    # read in txt file of counts
    counts = pd.read_csv(args.a, index_col="Gene", sep="\t")

    counts = processCounts(counts)

    counts = quantileNormalize(counts)

    counts = log2Transform(counts)

    counts_t = counts.T

    counts_t_new = counts_t.index.str.split(" ", n=1, expand=True)

    kayrotype = []
    for i in counts_t_new:
        kayrotype.append(i[1])

    gene_name = []
    for i in counts_t_new:
        gene_name.append(i[0])

    counts_t["Gene"] = gene_name
    counts_t["Kayrotype"] = kayrotype

    X = counts_t.iloc[:, 0:-2].values
    y = counts_t.iloc[:, -1].values


    colors = {'MLL': '#F6EE03',
              'Complex': '#F67403',
              '-7/(del(7)': '#C60006',
              'T(6:9)': '#0006C6',
              'Others': '#FCADFF',
              '-5/del(5)': '#5D0060'}

    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=3)

    Y_sklearn = sklearn_pca.fit_transform(X_std)

    data = []

    for name, col in zip(('MLL', 'Complex', '-7/(del(7)', 'T(6:9)', 'Others', '-5/del(5)'), colors.values()):


        trace = dict(
                type='scatter',
                x=Y_sklearn[y == name, 0],
                y=Y_sklearn[y == name, 2],
                    # z=Y_sklearn[y == name, 2],
                mode='markers',
                name=name,
                marker=dict(
                    color=col,
                    size=12,
                    line=dict(
                        color='rgba(0, 0, 0, 0.0)',
                        width=0.7),
                    opacity=0.7)
                )
        data.append(trace)

        layout = dict(
                xaxis=dict(title='PC1', showline=False),
                yaxis=dict(title='PC3', showline=False)
        )


    fig = go.Figure(data=data, layout=layout)

    fig.update_layout({'width': 600, 'height': 600,
                       'showlegend': True, 'hovermode': 'closest'})

    pca2d_div = plot(fig, include_plotlyjs=True, output_type='div')

    return pca2d_div



def scatterPlot2D_1_4():

    # read in txt file of counts
    counts = pd.read_csv(args.a, index_col="Gene", sep="\t")

    counts = processCounts(counts)

    counts = quantileNormalize(counts)

    counts = log2Transform(counts)

    counts_t = counts.T

    counts_t_new = counts_t.index.str.split(" ", n=1, expand=True)

    kayrotype = []
    for i in counts_t_new:
        kayrotype.append(i[1])

    gene_name = []
    for i in counts_t_new:
        gene_name.append(i[0])

    counts_t["Gene"] = gene_name
    counts_t["Kayrotype"] = kayrotype

    X = counts_t.iloc[:, 0:-2].values
    y = counts_t.iloc[:, -1].values


    colors = {'MLL': '#F6EE03',
              'Complex': '#F67403',
              '-7/(del(7)': '#C60006',
              'T(6:9)': '#0006C6',
              'Others': '#FCADFF',
              '-5/del(5)': '#5D0060'}

    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=4)

    Y_sklearn = sklearn_pca.fit_transform(X_std)

    data = []

    for name, col in zip(('MLL', 'Complex', '-7/(del(7)', 'T(6:9)', 'Others', '-5/del(5)'), colors.values()):


        trace = dict(
                type='scatter',
                x=Y_sklearn[y == name, 0],
                y=Y_sklearn[y == name, 3],
                    # z=Y_sklearn[y == name, 2],
                mode='markers',
                name=name,
                marker=dict(
                    color=col,
                    size=12,
                    line=dict(
                        color='rgba(0, 0, 0, 0.0)',
                        width=0.7),
                    opacity=0.7)
                )
        data.append(trace)

        layout = dict(
                xaxis=dict(title='PC1', showline=False),
                yaxis=dict(title='PC4', showline=False)
        )


    fig = go.Figure(data=data, layout=layout)

    fig.update_layout({'width': 600, 'height': 600,
                       'showlegend': True, 'hovermode': 'closest'})

    pca2d_div = plot(fig, include_plotlyjs=True, output_type='div')

    return pca2d_div




plot_3d = scatterPlot3D()
plot_2d = scatterPlot2D()
plot_2d_1_3 = scatterPlot2D_1_3()
plot_2d_1_4 = scatterPlot2D_1_4()


pca_3d.write(plot_3d)
pca_1_2_2d.write(plot_2d)
pca_1_3_2d.write(plot_2d_1_3)
pca_1_4_2d.write(plot_2d_1_4)


