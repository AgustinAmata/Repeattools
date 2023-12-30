import argparse
import sys
import traceback
from uuid import uuid1
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
    help_exclude = """If selected, it excludes from the violin plots
    unknown data and data belonging to other repetitive elements:
    Low_complexity, Satellite, Simple_repeat, LTR, DNA, rRNA, snRNA,
    ARTEFACT, tRNA-CR1, tRNA-Core-RTE, tRNA-Deu-RTE, tRNA-L1, tRNA-RTE
    and tRNA"""
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

    help_output_folder = """Output directory name for the violin
    and the box plots"""
    parser.add_argument("--output", "-o", type=Path,
                        help=help_output_folder, required=True)

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
    box_file = arguments.box
    groups_file = arguments.groups
    out_folder = arguments.output

    if not out_folder.exists():
        out_folder.mkdir()

    log_number = uuid1()
    log_fhand = open(out_folder / f"REPlotDivergence.{log_number}.log", "w")
    msg = f"Command used: {' '.join(sys.argv)}\n"
    if violin_dir:
        msg += f"Input directory for violin plots: {violin_dir.resolve()}\n"
    if box_file:
        msg += f"Input file for box plots: {box_file.resolve()}\n"
    msg += f"Output folder: {out_folder.resolve()}\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()

    out_violin = out_folder / f"Violin_plots_{log_number}.png"
    out_box = out_folder / f"Box_plots_{log_number}.png"

    try:
        if violin_dir:
            print(f"{'-'*10} Generating violin plots for divergence {'-'*10}")
            if not names_file:
                msg = f"File for the names of the species analyzed by RECollector was not provided\n"
                print(msg)
                log_fhand.write(msg)
                log_fhand.flush()
                log_fhand.close()
                sys.exit()

            with open(names_file) as names:
                filehand_species = read_names_file(names)
                analyzed_species = list(filehand_species.values())
                print("Read names of species file for violin plots")
                msg = f"Groups file location: {names_file.resolve()}\n"
                print(msg)
                log_fhand.write(msg)
                log_fhand.flush()

            files_list = list(violin_dir.glob("*"))
            if exclude:
                new_list = []
                cols_to_exclude = ["Unknown", "Low_complexity",
                                   "Satellite", "Simple_repeat",
                                   "rRNA", "snRNA",
                                   "ARTEFACT", "tRNA-CR1",
                                   "tRNA-Core-RTE", "tRNA-Deu-RTE",
                                   "tRNA-L1", "tRNA-RTE", "tRNA"]
                for file in files_list:
                    for col in cols_to_exclude:
                        if col in file.name:
                            break
                    else:
                        new_list.append(file)
                files_list = new_list
                msg = f"Excluded files: {', '.join(cols_to_exclude)}\n"
                print(msg)
                log_fhand.write(msg)
                log_fhand.flush()

            get_divergence_violins(files_list, tree_fpath, analyzed_species, out_violin)
            print(f"{'-'*10} Generated violin plots for divergence {'-'*10}")
            msg = f"Violin plots created at: {out_violin.resolve()}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()

        if box_file:
            print(f"{'-'*10} Generating box plots for {box_file.name} {'-'*10}")
            if not groups_file:
                msg = f"File for the groups for the box plots was not provided\n"
                print(msg)
                log_fhand.write(msg)
                log_fhand.flush()
                log_fhand.close()
                sys.exit()

            with open(groups_file) as groups:
                species_and_groups = read_names_file(groups)
                print("Read groups file for box plots")
                msg = f"Groups file location: {groups_file.resolve()}\n"
                print(msg)
                log_fhand.write(msg)
                log_fhand.flush()

            get_divergence_boxplots(box_file, species_and_groups, out_box)
            print(f"{'-'*10} Generated box plots for {box_file.name} {'-'*10}")
            msg = f"Box plots created at: {out_box.resolve()}\n"
            print(msg)
            log_fhand.write(msg)
            log_fhand.flush()

        log_fhand.close()

    except Exception as e:
        msg = f"{'*'*10} An error occurred. See traceback below {'*'*10}\n"
        print(msg)
        log_fhand.write(msg)
        log_fhand.write(traceback.format_exc())
        log_fhand.close()
        raise

if __name__ == "__main__":
    main()