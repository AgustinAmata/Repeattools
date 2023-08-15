import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ete3 import Tree, NodeStyle
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .plot_eteTree import plot_tree

def get_count_matrix_heatmap(matrix_df, out_file, group_dict, dendro=False):
    """Generate heatmap from TE count matrix.
    
    Parameters
    ----------
    matrix_df : `pandas.DataFrame`
        DataFrame in wide form, where columns are the different
        subcategories of the chosen category, and rows are species.
    
    out_file : output file

    group_dict : dictionary
        Contains the groups defined by the user. Keys: species;
        values: group.
    
    dendro: bool, default: False
        If True, the generated heatmap will contain a dendrogram
        for the species (rows).
    """
    group_df = matrix_df.index.map(group_dict)
    groups = group_df.unique()
    group_pal = sns.color_palette("YlOrBr", len(groups))
    group_lut = dict(zip(groups, group_pal))
    group_colors = group_df.map(group_lut)

    #Heatmap creation with previous standardization
    cm_heat = sns.clustermap(matrix_df, standard_scale=1,
                             cmap="mako", col_cluster=False,
                             row_colors=group_colors,
                             row_cluster=dendro)
    
    handles = [Patch(facecolor=group_lut[name]) for name in group_lut]
    plt.legend(handles, group_lut, title='Species',
               bbox_to_anchor=(0.5, 1),
               bbox_transform=plt.gcf().transFigure,
               loc='upper center')
    
    cm_heat.savefig(out_file, dpi=300)

def get_count_matrix_pca(matrix_df, out_file, group_dict, show_names=False):
    """Generate PCA plot from TE count matrix.

    Parameters
    ----------
    matrix_df : `pandas.DataFrame`
        DataFrame in wide form, where columns are the different
        subcategories of the chosen category, and rows are species.

    out_file : output file

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
    fig, ax = plt.subplots()
    sns.scatterplot(data=pca_df, x=pca_df.PC1, y=pca_df.PC2,
                    hue="group", ax=ax)
    ax.legend(loc= "center left", title="Groups",
              bbox_to_anchor=(1.1,0.5))
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

def get_divergence_violins(div_df, tree_fpath, out_file):
    """Generate violin plots given a long-form DataFrame.

    For a DataFrame containing divergence data from RECollector,
    generate violin plots for each species and for each subcategory
    in alphabetical order or along phylogenetic data.

    Parameters
    ----------
    div_df : pandas DataFrame
        Long-form DataFrame consisting of at least three columns:
        each row consists of a repeat, in which the first column
        contains the name of the species, the second column contains
        the subcategory of the selected category in RECollector
        (whose name is the overall name of the column) and the third
        column contains its divergence value (given as a float);
        further columns will be ignored.

    tree_fpath : path to a Newick tree file
        If provided, violin plots will be ordered according to the data
        given by the tree (which will also appear in the final figure).

    out_file : output file
    """
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

        fig.savefig(out_file, dpi=300)

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

        fig.savefig(out_file, dpi=300)
