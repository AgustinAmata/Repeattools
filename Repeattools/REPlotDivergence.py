import argparse
from pathlib import Path

from src.generate_plots import get_divergence_violins, get_divergence_boxplots
from src.utils import read_names_file

def argument_parser():
    desc = """Generates violin plots and box plots 
    based on RECollector's species divergence output.
    Violin plots are generated for all the species in the
    file and for each category, all in the same file (by default,
    species are ordered alphabetically, but an additional
    newick tree file can be added for specific order).
    On the other hand, box plots are generated for a specific
    category introduced by the user. For the box plots, a tab-separated
    file indicating in each line the species and the group it belongs
    to (for example, a taxonomical clade), so that different species
    can be grouped together."""
    parser = argparse.ArgumentParser(description=desc)

    help_divergence_violin= """Folder of the divergence
    file(s) of RECollector for the construction of violin
    plots"""
    parser.add_argument("--violin", "-v", type=Path, default=False,
                        help=help_divergence_violin, required=False)
    help_input_names_file = """Text file containing the names of all 
    the species analyzed by RECollector, that is, the same file
    required for RECollector."""
    parser.add_argument("--names", "-n", type=Path,
                        help=help_input_names_file, required=False)
    help_exclude = """If selected, it excludes from violin plots
    unknown data and data belonging to repetitive elements
    of low complexity, apart from rRNA"""
    parser.add_argument("--exclude", "-e", help=help_exclude,
                        action="store_true", default=False,
                        required=False)
    help_div_tree = """Optional tree file in newick format for the violin
    plot"""
    parser.add_argument("--tree", "-t", type=str, default=False,
                        help=help_div_tree, required=False)

    help_divergence_box = """RECollector divergence file
    for plotting box plots"""
    parser.add_argument("--box", "-b", type=Path, default=False,
                        help=help_divergence_box, required=False)
    help_box_group_file = """Tab-separated file with the species and the
    groups they belong to"""
    parser.add_argument("--groups", "-g", type=Path,
                        help=help_box_group_file, required=False)

    help_violin_output_name = """Output file name for violin plots"""
    parser.add_argument("--outvio", "-V", type=Path, 
                        help=help_violin_output_name, required=False)

    help_violin_output_name = """Output file name for box plots"""
    parser.add_argument("--outbox", "-B", type=Path,
                        help=help_violin_output_name, required=False)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    violin_dir = arguments.violin
    names_file = arguments.names
    exclude = arguments.exclude
    tree_fpath = arguments.tree
    outvio_fpath = arguments.outvio
    box_file = arguments.box
    groups_file = arguments.groups
    outbox_fpath = arguments.outbox

    if violin_dir:
        print(f"{'-'*10} Generating violin plots for divergence {'-'*10}")
        with open(names_file) as names:
            filehand_species = read_names_file(names)
            analyzed_species = list(filehand_species.values())
            print("Read names of species file for violin plots")

        files_list = list(violin_dir.glob("*"))
        if exclude:
            new_list = []
            for file in files_list:
                for col in ["Unknown", "Low_complexity", "Satellite", "Simple_repeat", "LTR", "DNA", "rRNA"]:
                    if col in file.name:
                        break
                else:
                    new_list.append(file)
            files_list = new_list

        get_divergence_violins(files_list, tree_fpath, analyzed_species, outvio_fpath)
        print(f"{'-'*10} Generated violin plots for divergence {'-'*10}")

    if box_file:
        print(f"{'-'*10} Generating box plots for {box_file.name} {'-'*10}")
        with open(groups_file) as groups:
            species_and_groups = read_names_file(groups)
            print("Read groups file for box plots")

        get_divergence_boxplots(box_file, species_and_groups, outbox_fpath)
        print(f"{'-'*10} Generated box plots for {box_file.name} {'-'*10}")

if __name__ == "__main__":
    main()