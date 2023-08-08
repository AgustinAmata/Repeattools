import builtins as __builtin__
from itertools import chain
from math import floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from ete3 import Tree
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .plot_eteTree import plot_tree

#Generate heatmap from TE count matrix
def get_count_matrix_heatmap(matrix_df, out_folder):
    #Filter out unknown data
    excluded = ["Unknown", "none", "{'none':'none'}"]
    known_df = matrix_df.drop(columns=excluded, errors="ignore")

    #Heatmap creation with previous standardization
    cm_heat = sns.clustermap(known_df, standard_scale=1,
                             cmap="mako", col_cluster=False)
    
    out_file = out_folder.joinpath("heatmap.png")
    cm_heat.savefig(out_file)

#Generate PCA from TE count matrix
def get_count_matrix_pca(matrix_df, out_folder, show_names, group_dict):
    #Filter out unknown data
    excluded = ["Unknown", "none", "{'none':'none'}"]
    known_df = matrix_df.drop(columns=excluded, errors="ignore")

    #Standardization of data and 2-component PCA creation
    scaled_matrix = StandardScaler().fit_transform(known_df)

    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(scaled_matrix)

    #Gather PCA data to plot
    pca_df = pd.DataFrame(pca_data, index=known_df.index,
                          columns=['PC1', 'PC2'])
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)

    #Assign each species to its group
    pca_df["group"] = pca_df.index.map(group_dict)

    #Plotting the results
    fig, ax = plt.subplots()
    sns.scatterplot(data=pca_df, x=pca_df.PC1, y=pca_df.PC2,
                    hue="group", ax=ax)
    ax.set_xlabel(f"PC1 ({per_var[0]}%)")
    ax.set_ylabel(f"PC2 ({per_var[1]}%)")

    #Show the species name of each point
    if show_names:
        for species in pca_df.index:
            ax.annotate(species, (pca_df.PC1.loc[species], pca_df.PC2.loc[species]))

    out_file = out_folder.joinpath("pca.png")
    fig.savefig(out_file)

#Generate violin plots for each species and for each category
#in alphabetical order or along phylogenetic data
def get_divergence_violins(div_df, tree_fpath, out_folder):
    #Filter out unknown data
    excluded = ["Unknown", "none", "{'none':'none'}"]
    div_df = div_df.drop(columns=excluded, errors="ignore")

    #No tree file provided
    if not tree_fpath:
        #Get categories
        cat_name = div_df.columns[1]
        cats = div_df[cat_name].unique()

        #Creating main figure and subplots
        width, height = plt.rcParams.get("figure.figsize")
        fig, axs = plt.subplots(1, len(cats),
                                figsize=(len(cats)*2, height*2),
                                sharey=True, sharex=True,
                                constrained_layout=True)
        #Create alphabetically ordered list
        n_species = sorted(div_df.species.unique())

        first_axes = True
        for cat, ax in zip(cats, axs):
            #Check if some species is not in the category
            cat_df = div_df[div_df[cat_name] == cat]
            sp_in_df = list(cat_df["species"].unique())
            for species in n_species:
                #Add blank data if not present
                if species not in sp_in_df:
                    empty_df = pd.DataFrame({"species": [species],
                                             cat_name: cat, "per div": 0.0})
                    cat_df = pd.concat([cat_df, empty_df])

            #Generate violin plot for the category,
            #limit the extent of the violin within the range
            #of the observed data
            sns.violinplot(data=cat_df, x="per div", y="species",
                           ax=ax, cut=0, order=n_species)

            # Hide the right, left, and top spines
            ax.set_title(cat)
            ax.set_ylabel("")
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)

            # Only show ticks on the left and bottom spines
            ax.xaxis.set_ticks_position("bottom")
            if first_axes:
                first_axes = False
            else:
                ax.yaxis.set_ticks_position("none")

        out_file = out_folder.joinpath("violins_alphabetical.png")
        fig.savefig(out_file)

    #Tree file provided
    else:
        #Read the Newick tree
        tree = Tree(tree_fpath)
        
        #Get categories
        cat_name = div_df.columns[1]
        cats = div_df[cat_name].unique()

        #Creating main figure and GridSpec
        width, height = plt.rcParams.get("figure.figsize")
        fig = plt.figure(figsize=(len(cats)*2+1, height*2),
                         constrained_layout=True)
        gs = GridSpec(1, len(cats)+2, figure=fig)

        #Tree takes two subplot spaces
        ax1 = fig.add_subplot(gs[0:2])
        coords = plot_tree(tree, axe=ax1, font_size=12)
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        xmax *= 1.2
        ymax *= 1.05

        #Get order of the species in the tree
        n_species = tree.get_leaf_names()
        #Create first Axes for violins so that they can share
        #a common x axis; all Axes are aligned with the tree
        ax2 = fig.add_subplot(gs[2], sharey=ax1)
        for cat, i in zip(cats, range(2, len(cats)+2)):
            #First Axes for the violins
            if i == 2:
                ax = ax2
            #Following Axes
            else:
                ax = fig.add_subplot(gs[i], sharey=ax1, sharex=ax2)
            
            #Check if some species is not in the category
            cat_df = div_df[div_df[cat_name] == cat]
            sp_in_df = list(cat_df["species"].unique())
            for species in n_species:
                #Add blank data if not present
                if species not in sp_in_df:
                    empty_df = pd.DataFrame({"species": [species],
                                             cat_name: cat, "per div": 0.0})
                    cat_df = pd.concat([cat_df, empty_df])

            #Generate violin plot for the category,
            #limit the extent of the violin within the range
            #of the observed data
            sns.violinplot(data=cat_df, x="per div", y="species",
                           ax=ax, cut=0, order=n_species[::-1])

            # Hide the right, left, and top spines
            ax.set_title(cat)
            ax.set_ylabel("")
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)

            # Only show ticks on the left and bottom spines
            ax.xaxis.set_ticks_position("bottom")
            ax.set_yticks([])

        _ = ax1.set_xlim(xmin, xmax)
        _ = ax1.set_ylim(ymin, ymax)

        out_file = out_folder.joinpath("violins_tree.png")
        fig.savefig(out_file)
