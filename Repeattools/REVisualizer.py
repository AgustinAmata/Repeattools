import argparse
import builtins as __builtin__
from itertools import chain
from math import floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pathlib import Path

from ete3 import Tree
from matplotlib.collections import LineCollection
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from src.generate_plots import (get_count_matrix_heatmap,
                                get_count_matrix_pca,
                                get_divergence_violins)
from src.utils import read_names_file, get_large_dfs
from src.plot_eteTree import plot_tree

def argument_parser():
    desc = """Generates plots based on RECollector's outputs.
    In the case of the TE count matrix, it generates a PCA and
    a heatmap whereas in the case of the species divergence file
    it generates a violin plot for all the species in the file
    (by default, the ridgeline plot is ordered alphabetically, but 
    an additional newick tree file can be added for specific order)"""
    parser = argparse.ArgumentParser(description=desc)

    help_matrix_input = """Input for the TE count matrix of 
    RECollector"""
    parser.add_argument("--matrix", "-m", type=Path, 
                        help=help_matrix_input, required=True)
    help_show_names = """If selected, points in the PCA plot will
    appear with their species name"""
    parser.add_argument("--names", "-n", help=help_show_names,
                        action="store_true", default=False,
                        required=False)
    help_group_file = """Include a file which includes in each
    line the name of the species and the selected taxon, separated
    by a tab"""
    parser.add_argument("--gfile", "-g", type=Path,
                        help=help_group_file, required=True)

    help_divergence_input = """Input for the species divergence
    file of RECollector"""
    parser.add_argument("--div", "-d", type=Path,
                        help=help_divergence_input, required=True)
    help_div_tree = """Optional tree file in newick format"""
    parser.add_argument("--tree", "-t", type=str, default=False,
                        help=help_div_tree, required=False)

    help_output_folder = """Output folder path"""
    parser.add_argument("--output", "-o", type=Path,
                        help=help_output_folder, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    matrix_fpath = arguments.matrix
    show_names = arguments.names
    group_fpath = arguments.gfile
    div_fpath = arguments.div
    tree_fpath = arguments.tree
    out_fpath = arguments.output

    with open(group_fpath) as gfile:
        group_dict = read_names_file(gfile)
        print("Read groups file")

    with open(matrix_fpath) as matrix:
        matrix_df = get_large_dfs(matrix, transpose=True)
        print("Read TE count matrix file")
        get_count_matrix_heatmap(matrix_df, out_fpath)
        print("Generated heatmap")
        get_count_matrix_pca(matrix_df, out_fpath,
                                 show_names, group_dict)
        print("Generated PCA")

    with open(div_fpath) as diver:
        div_df = get_large_dfs(diver)
        print("Read species divergence file")
        get_divergence_violins(div_df, tree_fpath, out_fpath)
        print("Generated violin plots for divergence")

if __name__ == "__main__":
    main()