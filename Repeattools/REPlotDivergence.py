import argparse
from pathlib import Path

import pandas as pd

from src.generate_plots import get_divergence_violins
from src.utils import get_large_dfs
from src.utils import read_names_file

def argument_parser():
    desc = """Generates plots based on RECollector's species
    divergence output. It generates violin plots for all the
    species in the file and for each category (by default,
    species are ordered alphabetically, but an additional
    newick tree file can be added for specific order)"""
    parser = argparse.ArgumentParser(description=desc)

    help_divergence_input = """Folder of the divergence
    file(s) of RECollector"""
    parser.add_argument("--input", "-i", type=Path,
                        help=help_divergence_input, required=True)
    help_input_names_file = """Text file containing the names of all 
    the species analyzed by RECollector, that is, the same file
    required for RECollector."""
    parser.add_argument("--names", "-n", type=Path,
                        help=help_input_names_file, required=True)
    help_exclude_unknown = """If selected, it excludes unknown data"""
    parser.add_argument("--exclude", "-e", help=help_exclude_unknown,
                        action="store_true", default=False,
                        required=False)
    help_div_tree = """Optional tree file in newick format"""
    parser.add_argument("--tree", "-t", type=str, default=False,
                        help=help_div_tree, required=False)

    help_output_name = """Output file name"""
    parser.add_argument("--output", "-o", type=Path,
                        help=help_output_name, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    div_dir = arguments.input
    names_file = arguments.names
    exclude = arguments.exclude
    tree_fpath = arguments.tree
    out_fpath = arguments.output

    with open(names_file) as names:
        filehand_species = read_names_file(names)
        analyzed_species = list(filehand_species.values())
        print("Read names of species file")
    
    files_list = list(div_dir.glob("*"))
    if exclude:
        files_list = list(filter(lambda x: "Unknown" not in x.name, files_list))

    print(f"{'-'*10} Generating violin plots for divergence {'-'*10}")
    get_divergence_violins(files_list, tree_fpath, analyzed_species, out_fpath)
    print(f"{'-'*10} Generated violin plots for divergence {'-'*10}")

if __name__ == "__main__":
    main()