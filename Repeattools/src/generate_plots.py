import builtins as __builtin__
from itertools import chain
from math import floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from ete3 import Tree
from matplotlib.collections import LineCollection
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .plot_eteTree import plot_tree

def get_count_matrix_heatmap(matrix_df, out_folder):
    excluded = ["Unknown", "none", "{'none':'none'}"]
    known_df = matrix_df.drop(columns=excluded, errors="ignore")
    cm_heat = sns.clustermap(known_df, standard_scale=1,
                             cmap="mako", col_cluster=False)
    out_file = out_folder.joinpath("heatmap.png")
    cm_heat.savefig(out_file)

def get_count_matrix_pca(matrix_df, out_folder, show_names, group_dict):
    excluded = ["Unknown", "none", "{'none':'none'}"]
    known_df = matrix_df.drop(columns=excluded, errors="ignore")
    scaled_matrix = StandardScaler().fit_transform(known_df)
    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(scaled_matrix)
    pca_df = pd.DataFrame(pca_data, index=known_df.index,
                          columns=['PC1', 'PC2'])
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)

    pca_df["group"] = pca_df.index.map(group_dict)

    fig, ax = plt.subplots()
    sns.scatterplot(data=pca_df, x=pca_df.PC1, y=pca_df.PC2,
                    hue="group", ax=ax)
    ax.set_xlabel(f"PC1 ({per_var[0]}%)")
    ax.set_ylabel(f"PC2 ({per_var[1]}%)")

    if show_names:
        for species in pca_df.index:
            ax.annotate(species, (pca_df.PC1.loc[species], pca_df.PC2.loc[species]))

    out_file = out_folder.joinpath("pca.png")
    fig.savefig(out_file)

def get_divergence_violins(div_df, tree_fpath, out_folder):
    excluded = ["Unknown", "none", "{'none':'none'}"]
    div_df = div_df.drop(columns=excluded, errors="ignore")

    if not tree_fpath:
        div_df = div_df.sort_values(by="species")
        fig, ax = plt.subplots()

        out_file = out_folder.joinpath("violins_alphabetical.png")
        fig.savefig(out_file)

    else:
        tree = Tree(tree_fpath)
        cat_name = div_df.columns[1]
        cats = div_df[cat_name].unique()
        width, height = plt.rcParams.get('figure.figsize')
        fig, axs = plt.subplots(1, len(cats)+1,
                                figsize=(len(cats)*2, height*2),
                                sharey=True, constrained_layout=True)
        ax1 = axs[0]
        coords = plot_tree(tree, axe=ax1)
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        xmax *= 1.2

        n_species = tree.get_leaf_names()

        for cat, ax in zip(cats, axs[1:]):
            cat_df = div_df[div_df[cat_name] == cat]
            sp_in_df = list(cat_df["species"].unique())
            for species in n_species:
                if species not in sp_in_df:
                    empty_df = pd.DataFrame({"species": [species],
                                             cat_name: cat, "per div": 0.0})
                    cat_df = pd.concat([cat_df, empty_df])

            sns.violinplot(data=cat_df, x="per div", y="species",
                           ax=ax, cut=0, sharey=ax1, order=n_species[::-1])

            # Hide the right and top spines
            ax.set_title(cat)
            ax.set_ylabel("")
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)

            # Only show ticks on the left and bottom spines
            ax.xaxis.set_ticks_position('bottom')
            ax.set_yticks([])

        _ = ax1.set_xlim(xmin, xmax)
        _ = ax1.set_ylim(ymin, ymax)

        out_file = out_folder.joinpath("violins_tree.png")
        fig.savefig(out_file)
