import argparse
from pathlib import Path

from src.generate_plots import (get_count_matrix_heatmap,
                                get_count_matrix_pca)
from src.utils import read_names_file, get_large_dfs

def argument_parser():
    desc = """Generates plots based on RECollector's TE count
    matrix output. It generates a PCA and a heatmap. A file
    detailing the group of each species must be provided
    (file in which each row corresponds to a species:
    the name of the species will appear first, followed by
    the name of the group it belongs, separated by a tab)."""
    parser = argparse.ArgumentParser(description=desc)

    help_matrix_input = """Input for the TE count matrix of 
    RECollector"""
    parser.add_argument("--input", "-i", type=Path, 
                        help=help_matrix_input, required=True)
    help_exclude_unknown = """If selected, it excludes unknown data"""
    parser.add_argument("--exclude", "-e", help=help_exclude_unknown,
                        action="store_true", default=False,
                        required=False)
    help_show_dendrogram = """If selected, heatmap will contain
    a dendrogram for the species"""
    parser.add_argument("--dendro", "-d", help=help_show_dendrogram,
                        action="store_true", default=False,
                        required=False)
    help_show_names = """If selected, points in the PCA plot will
    include their species name"""
    parser.add_argument("--names", "-n", help=help_show_names,
                        action="store_true", default=False,
                        required=False)
    help_group_file = """Include a file which includes in each
    line the name of the species and the selected group, separated
    by a tab"""
    parser.add_argument("--gfile", "-g", type=Path,
                        help=help_group_file, required=True)

    help_output_heatmap = """Output file name for the heatmap"""
    parser.add_argument("--heatmap", "-H", type=Path,
                        help=help_output_heatmap, required=True)
    help_output_pca = """Output file name for the PCA"""
    parser.add_argument("--pca", "-p", type=Path,
                        help=help_output_pca, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    matrix_fpath = arguments.matrix
    exclude = arguments.exclude
    dendro = arguments.dendro
    show_names = arguments.names
    group_fpath = arguments.gfile
    out_heatmap = arguments.heatmap
    out_pca = arguments.pca

    with open(group_fpath) as gfile:
        group_dict = read_names_file(gfile)
        print("Read groups file")

    with open(matrix_fpath) as matrix:
        matrix_df = get_large_dfs(matrix, exclude, transpose=True)
        print("Read TE count matrix file")
        get_count_matrix_heatmap(matrix_df, out_heatmap,
                                 group_dict, dendro=dendro)
        print("Generated heatmap")
        get_count_matrix_pca(matrix_df, out_pca,
                             group_dict, show_names=show_names)
        print("Generated PCA")

if __name__ == "__main__":
    main()