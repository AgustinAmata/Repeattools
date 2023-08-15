import argparse
from pathlib import Path


from src.create_matrix import (create_df_from_parsed_input,
                               create_te_count_matrix,
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
    TESorter files must come from RM. Several options are provided
    to filter the data from these files (some of them require input
    files). Another file with the names of the species is required."""
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
    by a tab. Example:Peame105,Persea_americana"""
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
    (class, superfamily, tes order, tes superfamily, clade, domains)"""
    depth_choices = ["superfamily","class", "tes order",
                     "tes superfamily", "clade", "domains"]
    parser.add_argument("--depth", help=help_matrix_depth,
                        choices=depth_choices,
                        default="superfamily", required=False)

    help_output = """Output folder path. Generated files will
    be in .csv format"""
    parser.add_argument("--output", "-o", help=help_output,
                        required=True)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def main():
    arguments = get_options()

    root_dir = arguments.input
    rm_dir = root_dir / "RepeatMasker"
    te_dir = root_dir / "TESorter"
    rm_files = rm_dir.glob("*.out")
    te_files = te_dir.glob("*.tsv")

    names_file = arguments.names
    with open(names_file) as names:
        filehand_species = read_names_file(names)
        print("Read names of species file")

    read_inputs = {}

    for rm_file in rm_files:
        with open(rm_file) as rm_fhand:
            fhand_repeats = read_repeatmasker_out(rm_fhand)
            for filehand in filehand_species.keys():
                if filehand in rm_file.name:
                    sp_name = filehand_species[filehand]
                    read_inputs[sp_name] = [fhand_repeats]
                    break
            print(f"Read {rm_file.name} -> {sp_name}")

    for te_file in te_files:
        with open(te_file) as te_fhand:
            fhand_repeats = read_tesorter_cls_tsv(te_fhand)
            for filehand in filehand_species.keys():
                if filehand in te_file.name:
                    sp_name = filehand_species[filehand]
                    read_inputs[sp_name].append(fhand_repeats)
                    break
            print(f"Read {te_file.name} -> {sp_name}")

    merged_inputs = {}
    for inp in read_inputs:
        rm_input = read_inputs[inp][0]
        te_input = read_inputs[inp][1]
        merged_inp = merge_inputs(rm_input, te_input)
        merged_inputs[inp] = merged_inp
        print(f"Merged inputs from {inp}")

    print("Creating dataframes")
    species_dfs = {}
    for species in merged_inputs:
        sp_rep_list = merged_inputs[species]
        sp_df = create_df_from_parsed_input(sp_rep_list)
        species_dfs[species] = sp_df
    print("Dataframes created")

    if arguments.length:
        print("Started filtering by length")
        length = arguments.length
        for species in species_dfs:
            sp_df = species_dfs[species]
            filtered_df = filter_df_by_length(sp_df, length)
            species_dfs.update({species:filtered_df})
            print(f"Filtered data from {species}")

        print("Finished filtering by length")

    if arguments.chrom:
        print("Started filtering by chromosomes")

        chroms_file = arguments.chrom
        with open(chroms_file) as chroms:
            chrs_to_filter = read_chroms_file(chroms)
            print("Read file detailing chromosomes to filter")

        for species in chrs_to_filter:
            sp_df = species_dfs[species]
            selected_chrs = chrs_to_filter[species]
            chr_mode = arguments.E
            filtered_df = filter_df_by_chromosomes(sp_df, selected_chrs, chr_mode)
            species_dfs.update({species:filtered_df})
            print(f"Filtered data from {species}")

        print("Finished filtering by chromosomes")

    if arguments.domains:
        print("Started filtering by domains")
        if arguments.D:
            with open(arguments.D) as doms_file:
                domains, clades, features_dict = read_doms_file(doms_file)
            print("Read data from additional file")

        else:
            domains = []
            clades = []
            features_dict = []

        for species in species_dfs:
            sp_df = species_dfs[species]
            filtered_df = filter_df_by_domain(sp_df, domains, clades, features_dict)
            species_dfs.update({species:filtered_df})
            print(f"Filtered data from {species}")

        print("Finished filtering by domains")

    if arguments.per:
        print("Started filtering by percentage")

        threshold = arguments.t
        perc_mode = arguments.m
        for species in species_dfs:
            sp_df = species_dfs[species]
            filtered_df = filter_df_by_percentages(sp_df, threshold, perc_mode)
            species_dfs.update({species:filtered_df})
            print(f"Filtered data from {species}")

        print("Finished filtering by percentage")

    species_counted_tes = []
    depth = arguments.depth
    for species in species_dfs:
        sp_df = species_dfs[species]
        counted_tes = count_tes(sp_df, species, depth)
        species_counted_tes.append(counted_tes)

    print("Creating TE count matrix")
    te_count_matrix = create_te_count_matrix(species_counted_tes)
    print("TE count matrix created")

    print("Creating species divergence data")
    long_df_div = convert_data_to_long_df_div(species_dfs, depth)
    print("Species divergence data created")

    out_folder = Path(arguments.output)
    if not out_folder.exists():
        out_folder.mkdir()

    c_matrix_fpath = out_folder.joinpath(f"{out_folder.name}_count_matrix.csv")
    div_csv_fpath = out_folder.joinpath(f"{out_folder.name}_divergence.csv")

    te_count_matrix.to_csv(c_matrix_fpath, index_label= depth)
    print("TE count matrix file generated")

    long_df_div.to_csv(div_csv_fpath)
    print("Species divergence data file generated")

if __name__ == "__main__":
    main()