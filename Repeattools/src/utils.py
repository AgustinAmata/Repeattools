from csv import DictReader
import pandas as pd

def convert_data_to_long_df_div(species_dfs, depth):
    long_df_div = pd.DataFrame()
    for species in species_dfs:
        sp_df = species_dfs[species]
        sp_df.insert(0, "species", species)
        selected_columns = ["species", depth, "per div"]
        long_df = sp_df.loc[:, selected_columns]
        long_df_div = pd.concat([long_df_div, long_df], axis= 0)
    long_df_div = long_df_div.reset_index(drop=True)
    return long_df_div

def get_large_dfs(file, transpose=False):
    df_chunk = pd.read_csv(file, header=0, index_col=0, chunksize=1000000)
    chunk_list = []
    for chunk in df_chunk:
        chunk_list.append(chunk)

    df_concat = pd.concat(chunk_list)
    if transpose:
        df_concat = df_concat.T

    return df_concat

#Reads the file for chromosome filtering
def read_chroms_file(chr_file):
    chrs_to_filter = {}
    for species_chrom in chr_file:
        species_chrom = species_chrom.strip("\n").split("\t")
        sp_name = species_chrom[0]
        chrs = species_chrom[1].split(",")
        chrs_to_filter[sp_name] = chrs
    return chrs_to_filter

#Reads the file for domain filtering
def read_doms_file(doms_file):
    dfile_data =  DictReader(doms_file, delimiter="\t")
    for row in dfile_data:
        domains = row["domains"]
        clades = row["clades"]
        features = row["features"]

        if domains:
            domains = domains.split(",")
        else:
            domains = []

        if clades:
            clades = clades.split(",")
        else:
            clades = []

        if features:
            features_list = features.split(",")
            features_dict = []
            for feat in features_list:
                feat = feat.split(":")
                fdom = feat[0]
                fclade = feat[1]
                dom_clade = {fdom:fclade}
                features_dict.append(dom_clade)
        else:
            features_dict = []

    return domains, clades, features_dict

#Reads the file for species name-file name association,
#also applicable for group name file in REVisualizer
def read_names_file(names_file):
    filehand_species = {}
    for pair in names_file:
        pair = pair.strip("\n").split("\t")
        filehand = pair[0]
        sp_name = pair[1]
        filehand_species[filehand] = sp_name
    return filehand_species