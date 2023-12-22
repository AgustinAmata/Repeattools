import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ete3 import Tree, NodeStyle
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .plot_eteTree import plot_tree
from .utils import get_large_dfs

def get_count_matrix_heatmap(matrix_df, out_file, group_dict, hsize, dendro=False):
    """Generate heatmap from TE count matrix.
    
    Parameters
    ----------
    matrix_df : `pandas.DataFrame`
        DataFrame in wide form, where columns are the different
        subcategories of the chosen category, and rows are species.
    
    out_file : output file path

    group_dict : dictionary
        Contains the groups defined by the user. Keys: species;
        values: group.

    hsize : tuple (width, height)
        Size of the figure in inches

    dendro: bool, default: False
        If True, the generated heatmap will contain a dendrogram
        for the species (rows).
    """
    #Create colors for the groups and the legend
    group_df = matrix_df.index.map(group_dict)
    groups = group_df.unique()
    group_pal = sns.color_palette("hls", len(groups))
    group_lut = dict(zip(groups, group_pal))
    group_colors = group_df.map(group_lut)

    #Heatmap creation with previous standardization
    standard_matrix = matrix_df.apply(lambda col: zscore(col) if col.std() != 0 else col, axis=0)
    cm_heat = sns.clustermap(standard_matrix,
                             cmap="magma_r", col_cluster=True,
                             row_colors=group_colors,
                             row_cluster=dendro, figsize=hsize,
                             cbar_kws={"orientation":"horizontal"}, vmin=0)
    cm_heat.ax_col_dendrogram.set_visible(False)
    dendro_box = cm_heat.ax_col_dendrogram.get_position()
    dendro_box.y0 = dendro_box.y0 + 0.02
    dendro_box.y1 = dendro_box.y0 + 0.02
    cm_heat.cax.set_position(dendro_box)

    #Add legend
    handles = [Patch(facecolor=group_lut[name]) for name in group_lut]
    fontsize_scaled = math.ceil(hsize[1])
    legend_box = cm_heat.ax_col_dendrogram.get_position()
    yper = 1.5/hsize[1]
    legend_box.y0 = dendro_box.y1 + yper
    legend_box.y1 = legend_box.y0 + 0.08
    plt.legend(handles, group_lut, title='Groups',
               bbox_to_anchor=legend_box,
               ncol=math.ceil(len(handles)/4),
               bbox_transform=plt.gcf().transFigure,
               loc='upper center', title_fontsize=fontsize_scaled,
               fontsize=fontsize_scaled)
    cm_heat.ax_heatmap.set_xlabel("")
    cm_heat.savefig(out_file, dpi=200)

def get_count_matrix_pca(matrix_df, out_file, group_dict, show_names=False):
    """Generate PCA plot from TE count matrix.

    Parameters
    ----------
    matrix_df : `pandas.DataFrame`
        DataFrame in wide form, where columns are the different
        subcategories of the chosen category, and rows are species.

    out_file : output file path

    group_dict : dictionary
        Contains the groups defined by the user. Keys: species;
        values: group.
    
    show_names : bool, default: False
        If True, the generated plot will show the name of each point.
    """
    #Standardization of data and 2-component PCA creation
    scaled_matrix = StandardScaler().fit_transform(matrix_df)

    pca = PCA(n_components=2)
    pca_data = pca.fit_transform(scaled_matrix)

    #Gather PCA data to plot
    pca_df = pd.DataFrame(pca_data, index=matrix_df.index,
                          columns=['PC1', 'PC2'])
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)

    #Assign each species to its group
    pca_df["group"] = pca_df.index.map(group_dict)

    #Plotting the results
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    sns.scatterplot(data=pca_df, x=pca_df.PC1, y=pca_df.PC2,
                    hue="group", ax=ax)
    ax.legend(loc="center left", title="Groups",
              fontsize=8,
              bbox_to_anchor=(1, 0.5))
    ax.set_xlabel(f"PC1 ({per_var[0]}%)")
    ax.set_ylabel(f"PC2 ({per_var[1]}%)")

    ax.set_axisbelow(True)
    ax.grid(alpha=0.7)
    ax.axvline(linestyle=(5, (10, 3)), linewidth=1, color="gray", zorder=0.6)
    ax.axhline(linestyle=(5, (10, 3)), linewidth=1, color="gray", zorder=0.6)

    [ax.spines[side].set_visible(False) for side in ax.spines]
    #Show the species name of each point
    if show_names:
        for species in pca_df.index:
            x_pos = pca_df.PC1.loc[species]
            y_pos = pca_df.PC2.loc[species]
            trans = ax.transData.transform((x_pos, y_pos))
            tx, ty = ax.transAxes.inverted().transform(trans)
            if tx > 0.9:
                ha = "right"
            else:
                ha = "left"
            if ty > 0.9:
                va = "top"
            else:
                va = "bottom"
            ax.annotate(species, (x_pos, y_pos), ha=ha, va=va,
                        fontsize=7)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300)

def get_divergence_boxplots(div_file, species_and_groups, out_fpath):
    """Generate a box plot from a given category from RECollector divergence data.

    Parameters
    ----------
    div_file : path to the RECollector divergence file

    species_and_groups : dictionary
        Contains the groups defined by the user. Keys: species;
        values: group.

    out_fpath : output file path
    """
    sp_per_group = {}
    for k, v in species_and_groups.items():
        if v in sp_per_group:
            sp_per_group[v].append(k)
        else:
            sp_per_group[v] = [k]
    ordered_sp = [sp for val in sp_per_group.values() for sp in val]

    with open(div_file) as file:
        div_df = get_large_dfs(file)
        cat_name = div_df.columns[1]
        cat = list(div_df[cat_name].unique())[0]
        div_df["group"] = div_df["species"].apply(lambda x: species_and_groups[x]).astype("category")
        print(f"Read data for {cat} divergence")
    
    #Check if some species is not in the dataframe
    sp_in_df = list(div_df["species"].unique())    
    for species in ordered_sp:
        #Add blank data if not present
        if species not in sp_in_df:
            empty_df = pd.DataFrame({"species": [species],
                                        cat_name: cat, "per div": 0.0})
            div_df = pd.concat([div_df, empty_df])

    #Modify plot aspect given the number of species
    data_length = len(ordered_sp)
    if data_length <= 10:
        plot_aspect = 2
    elif data_length > 10 and data_length <= 25:
        plot_aspect = 3
    elif data_length > 25 and data_length <= 50:
        plot_aspect = 3.5
    elif data_length > 50 and data_length <= 75:
        plot_aspect = 4
    else:
        plot_aspect = 5

    b_plot = sns.catplot(data=div_df, x="species",
                         y="per div", hue="group",
                         kind="box", aspect=plot_aspect, order=ordered_sp,
                         hue_order=list(sp_per_group.keys()),
                         dodge=False, sharex=True)
    b_plot.set_xticklabels(rotation=90)
    sns.despine(bottom=True)
    b_plot.set(title=f"Divergence data for {cat}", 
               ylabel="% Divergence", xlabel="",
               axisbelow=True)
    b_plot.legend.set(title="Groups")
    b_plot.ax.yaxis.grid()
    b_plot.savefig(out_fpath, dpi=200)

def get_divergence_violins(files_list, tree_fpath, analyzed_species, out_file):
    """Generate violin plots given a long-form DataFrame.

    For a DataFrame containing divergence data from RECollector,
    generate violin plots for each species and for each subcategory
    in alphabetical order or along phylogenetic data.

    Parameters
    ----------
    files_list : list of paths
        List composed of the paths to the divergence files from RECollector.

    tree_fpath : path to a Newick tree file
        If provided, violin plots will be ordered according to the data
        given by the tree (which will also appear in the final figure).

    analyzed_species: list
        List containing the names of the species to analyze
    
    out_file : output file path
    """

    #No tree file provided
    if not tree_fpath:
        #Creating main figure and subplots
        width, height = plt.rcParams.get("figure.figsize")
        fig, axs = plt.subplots(1, len(files_list),
                                figsize=(len(files_list)*2, height*4),
                                sharey=True, sharex=True,
                                constrained_layout=True)
        #Create alphabetically ordered list
        n_species = sorted(analyzed_species)
        ax_count = 0
        first_axes = True
        #Create violin plots for each subcategory
        for i, file in enumerate(files_list):
            #Create DataFrame and check if some species is not it
            if i > 0:
                for n_ax in axs[:i]:        
                    if not n_ax.has_data():    
                        ax = n_ax
                        break
                else:
                    ax = axs[i]
            else:
                ax = axs[i]

            with open(file) as div_file:
                div_df = get_large_dfs(div_file)
                cat_name = div_df.columns[1]
                cat = list(div_df[cat_name].unique())[0]
                print(f"Read data for {cat} divergence")
            sp_in_df = list(div_df["species"].unique())
            if (len(sp_in_df)/len(n_species)) < 0.75:
                print("Not enough species to proceed with the plot")
                continue
            ax_count += 1
            for species in n_species:
                #Add blank data if not present
                if species not in sp_in_df:
                    empty_df = pd.DataFrame({"species": [species],
                                             cat_name: cat, "per div": 0.0})
                    div_df = pd.concat([div_df, empty_df])

            #Generate violin plot for the category,
            #limit the extent of the violin within the range
            #of the observed data
            sns.violinplot(data=div_df, x="per div", y="species",
                           ax=ax, cut=0, order=n_species)
            print(f"Generated violin plots for {cat} divergence")

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
        for ax in axs[ax_count:]:
            fig.delaxes(ax)
        # fig.set_size_inches(len(ax_count)*2, height*4)
        fig.savefig(out_file, bbox_inches="tight", dpi=200)

    #Tree file provided
    else:
        #Read the Newick tree
        tree = Tree(tree_fpath)
        nstyle = NodeStyle()
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_width"] = 2
        nstyle["size"] = 7
        for n in tree.traverse():
            n.set_style(nstyle)

        #Creating main figure and GridSpec
        width, height = plt.rcParams.get("figure.figsize")
        fig = plt.figure(figsize=(len(files_list)*2+1, height*2),
                         constrained_layout=True)
        gs = GridSpec(1, len(files_list)+2, figure=fig)

        #Tree takes two subplot spaces
        ax1 = fig.add_subplot(gs[0:2])
        coords = plot_tree(tree, axe=ax1, name_offset= 0.05, 
                           font_size=12)
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        xmax *= 1.2
        ymax *= 1.05

        #Get order of the species in the tree
        n_species = tree.get_leaf_names()
        #Create first Axes for violins so that they can share
        #a common x axis; all Axes are aligned with the tree
        ax2 = fig.add_subplot(gs[2], sharey=ax1)
        for file, i in zip(files_list, range(2, len(files_list)+2)):
            #First Axes for the violins
            if i == 2:
                ax = ax2
            #Following Axes
            else:
                ax = fig.add_subplot(gs[i], sharey=ax1, sharex=ax2)
            
            #Create DataFrame and check if some species is not it
            with open(file) as div_file:
                div_df = get_large_dfs(div_file)
                cat_name = div_df.columns[1]
                cat = list(div_df[cat_name].unique())[0]
                print(f"Read data for {cat} divergence")
            sp_in_df = list(div_df["species"].unique())
            for species in n_species:
                #Add blank data if not present
                if species not in sp_in_df:
                    empty_df = pd.DataFrame({"species": [species],
                                             cat_name: cat, "per div": 0.0})
                    div_df = pd.concat([div_df, empty_df])

            #Generate violin plot for the category,
            #limit the extent of the violin within the range
            #of the observed data
            sns.violinplot(data=div_df, x="per div", y="species",
                           ax=ax, cut=0, order=n_species[::-1])
            print(f"Generated violin plots for {cat} divergence")

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

        fig.savefig(out_file, dpi=200)
