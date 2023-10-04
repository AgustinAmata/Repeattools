import argparse
import gc
import sys
from pathlib import Path

from src.create_matrix import (create_te_count_matrix,
                               count_tes, filter_df_by_chromosomes,
                               filter_df_by_domain,
                               filter_df_by_length,
                               filter_df_by_percentages) 
from src.read_input import (merge_inputs, read_repeatmasker_out,
                            read_tesorter_cls_tsv)
from src.utils import (convert_data_to_long_df_div,
                       read_chroms_file, read_doms_file,
                       read_names_file)


def argument_parser():
    desc = """Create a TE count matrix and a divergence table 
    from several files of RepeatMasker (RM) and TESorter (TES);
    TES files must come from RM individual repeats. Several options
    are provided to filter the data from these files (some of them
    require input files). Another file with the names of the species
    is required."""
    parser = argparse.ArgumentParser(description=desc)
        
    help_input_dir = """Input directory containing the files
    from RM and TES. Directory must contain two subdirectories
    named RepeatMasker and TESorter"""
    parser.add_argument("--input", "-i", type=Path, help=help_input_dir,
                        required=True)

    help_input_names_file = """Text file containing the names of the species 
    for each RM/TES file. File must contain the common part for
    both RM and TES files, followed by the name of the species 
    (use underscores instead of whitespaces) and separated
    by a tab. Example: Peame105  Persea_americana"""
    parser.add_argument("--names", "-n", type=Path,
                        help=help_input_names_file, required=True)

    help_filter_len = """Elements with less than the specified length
    will be filtered out from the data"""
    parser.add_argument("--length", "-l", type=int,
                        help= help_filter_len, required=False)

    desc_filter_chrom = """If selected, data will be filtered to only
    contain (or exclude) the chromosomes of the species in the file.
    Requires a text file composed of the name of the species and the
    desired, comma-separated chromosomes (name and chrs must be
    separated with a tab). Example: Persea_americana    chr1,chr2,chr3"""
    chrom_group = parser.add_argument_group("Filter by chromosomes",
                                            description=desc_filter_chrom)
    help_filter_chrom = """Select to allow filtering by chromosomes.
    By default data will only contain the selected chromosomes
    for the target species"""
    chrom_group.add_argument("--chrom", "-c", type=Path,
                             help=help_filter_chrom, required=False)
    help_exclude_chrom = """Select to exclude the specified chromosomes
    from the data"""
    chrom_group.add_argument("-E", help=help_exclude_chrom, default=False,
                             action="store_true", required=False)

    desc_filter_domains = """By default, it filters data to include
    repeats that contain specified domains (given by TES).
    An additional file can be introduced to further filter
    by specific domains, clades or domain:clade features."""
    dom_group = parser.add_argument_group("Filter by domains",
                                          description=desc_filter_domains)
    help_filter_domains = """Select to allow to allow the default filtering.
    Also required for the other option"""
    dom_group.add_argument("--domains", "-d", help=help_filter_domains,
                        action="store_true", required=False)
    help_specify_filter = """Introduce additional file to further filter
    by specific domains, clades or domain:clade features. First line must be:
    'domains    clades  features' (separated by tabs).
    Second line will be contain the values to filter in each
    category, separated by tabs; values in the same category
    must be separated by commas"""
    dom_group.add_argument("-D", type=Path, help=help_specify_filter,
                           default=[], required=False)

    desc_filter_percentages = """By default, it will filter the data
    to only include repeats with a percentage of divergence lower than 20%"""
    perc_group = parser.add_argument_group("Filter by percentages",
                                          description=desc_filter_percentages)
    help_filter_perc = """Select to allow filtering by percentages.
    Also required for the other options"""
    perc_group.add_argument("--per", "-p", help=help_filter_perc,
                        action="store_true",required=False)
    help_filter_perc_threshold = """Select the percentage of divergence,
    deletions or insertions to use as a threshold. Percentage must be
    preceded by div/del/ins= (e.g. del=30.0)"""
    perc_group.add_argument("-t", help=help_filter_perc_threshold,
                            default="div=20.0", nargs=1, required=False)
    help_filter_perc_mode = """Filters data to only include
    repeats lower, higher or equal to the threshold"""
    perc_group.add_argument("-m", help=help_filter_perc_mode,
                            default="lower_than",
                            choices=["lower_than","higher_than","equal"],
                            nargs=1)

    help_matrix_depth = """Select the depth of the TE count matrix
    (class, superfamily, tes order, tes superfamily, clade).
    Default superfamily."""
    depth_choices = ["superfamily", "class", "tes order",
                     "tes superfamily", "clade"]
    parser.add_argument("--depth", help=help_matrix_depth,
                        choices=depth_choices,
                        default="superfamily", required=False)
    help_override = """When selected, and when selected depth for the TE
    count matrix is 'superfamily' or 'class', repeats with unknown values
    for these columns will be overriden by their correspondent
    'tes superfamily' and 'tes order', respectively (if there are values
    for these columns in a given repeat)."""
    parser.add_argument("--override", help=help_override,
                        action="store_true", required=False)
    help_output = """Output folder path. Generated files will
    be in .csv format"""
    parser.add_argument("--output", "-o", help=help_output,
                        type=Path, required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()
    depth = arguments.depth
    override = arguments.override

    root_dir = arguments.input
    rm_dir = root_dir / "RepeatMasker"
    te_dir = root_dir / "TESorter"

    out_folder = arguments.output
    if not out_folder.exists():
        out_folder.mkdir()
    div_folder = out_folder.joinpath(f"{depth}_divergence_files")
    if not div_folder.exists():
        div_folder.mkdir()

    names_file = arguments.names
    with open(names_file) as names:
        filehand_species = read_names_file(names)
        print("Read names of species file")

    if arguments.chrom:
        chroms_file = arguments.chrom
        chr_mode = arguments.E
        with open(chroms_file) as chroms:
            chrs_to_filter = read_chroms_file(chroms)
            print("Read file detailing chromosomes to filter")

    if arguments.domains:
        if arguments.D:
            with open(arguments.D) as doms_file:
                domains, clades, features_dict = read_doms_file(doms_file)
            print("Read data from additional file")

        else:
            domains = []
            clades = []
            features_dict = []
            print("Created conditions for domain filtering")

    if arguments.per:
        threshold = arguments.t
        perc_mode = arguments.m
        print("Gathered conditions for percentage filtering")

    species_counted_tes = []
    for filehand in filehand_species:
        species = filehand_species[filehand]
        print(f"{'-'*10} Collecting data for {species} {'-'*10}")
        rm_file = list(rm_dir.glob(f"{filehand}*.out"))
        te_file = list(te_dir.glob(f"{filehand}*.cls.tsv"))

        if len(rm_file) != 1:
            print("RepeatMasker file was not found/file name in names file was not clear enough")
            print("Proceding with next species")
            continue

        if len(te_file) != 1:
            print("RepeatMasker file was not found/file name in names file was not clear enough")
            print("Proceding with next species")
            continue

        with open(rm_file[0]) as rm_fhand:
            print(f"Reading {rm_file[0].name}")
            rm_repeats = read_repeatmasker_out(rm_fhand)
            print(f"Read {rm_file[0].name}")
        with open(te_file[0]) as te_fhand:
            print(f"Reading {te_file[0].name}")    
            te_repeats = read_tesorter_cls_tsv(te_fhand)
            print(f"Read {te_file[0].name}")

        species_df = merge_inputs(rm_repeats, te_repeats)
        print("Merged input files into a dataframe")

        del rm_repeats, te_repeats

        if arguments.length:
            print("Started filtering by length")
            length = arguments.length
            species_df = filter_df_by_length(species_df, length)
            print("Finished filtering by length")

        if arguments.chrom:
            print("Started filtering by chromosomes")
            selected_chrs = chrs_to_filter[species]
            species_df = filter_df_by_chromosomes(species_df,
                                                  selected_chrs,
                                                  chr_mode)
            print("Finished filtering by chromosomes")

        if arguments.domains:
            print("Started filtering by domains")
            species_df = filter_df_by_domain(species_df,
                                             domains, clades,
                                             features_dict)
            print("Finished filtering by domains")

        if arguments.per:
            print("Started filtering by percentage")
            species_df = filter_df_by_percentages(species_df,
                                                  threshold,
                                                  perc_mode)
            print("Finished filtering by percentage")

        print(f"Counting TEs for {depth}")
        counted_tes = count_tes(species_df, species, depth, override)
        species_counted_tes.append(counted_tes)
        print(f"Counted TEs for {depth}")
        
        print(f"{'*'*5} Creating {depth} divergence data file(s) {'*'*5}") 
        depth_cats = species_df[depth].unique()   
        for cat in depth_cats:
            cat_df = species_df.loc[species_df[depth] == cat]
            long_df_div = convert_data_to_long_df_div(cat_df, species, depth)
            div_csv_fpath = div_folder.joinpath(f"{cat}_divergence.csv")
            if not div_csv_fpath.exists():
                long_df_div.to_csv(div_csv_fpath, index=False,
                                    chunksize=100000)
                print(f"{cat.capitalize()} divergence data file created")

            else:
                long_df_div.to_csv(div_csv_fpath, mode="a",
                                    index=False, header=False,
                                    chunksize=100000)
                print(f"{cat.capitalize()} divergence data file updated")

        del species_df, long_df_div
        gc.collect()
    
    print(f"{'-'*10} Performed operations for all accepted species {'-'*10}")
    print("Creating TE count matrix")
    te_count_matrix = create_te_count_matrix(species_counted_tes)
    print("TE count matrix created")

    c_matrix_fpath = out_folder.joinpath(f"{out_folder.name}_count_matrix.csv")

    te_count_matrix.to_csv(c_matrix_fpath, index_label=depth)
    print("TE count matrix file generated")

if __name__ == "__main__":
    main()