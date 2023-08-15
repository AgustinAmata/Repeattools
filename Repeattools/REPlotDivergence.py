import argparse
from pathlib import Path

from src.generate_plots import get_divergence_violins
from src.utils import get_large_dfs

def argument_parser():
    desc = """Generates plots based on RECollector's species
    divergence output. It generates violin plots for all the
    species in the file and for each category (by default,
    species are ordered alphabetically, but an additional
    newick tree file can be added for specific order)"""
    parser = argparse.ArgumentParser(description=desc)

    help_divergence_input = """Input for the species divergence
    file of RECollector"""
    parser.add_argument("--input", "-i", type=Path,
                        help=help_divergence_input, required=True)
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
    div_fpath = arguments.div
    exclude = arguments.exclude
    tree_fpath = arguments.tree
    out_fpath = arguments.output

    with open(div_fpath) as diver:
        div_df = get_large_dfs(diver, exclude)
        print("Read species divergence file")
        get_divergence_violins(div_df, tree_fpath, out_fpath)
        print("Generated violin plots for divergence")

if __name__ == "__main__":
    main()