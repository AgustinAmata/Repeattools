from csv import DictReader
import pandas as pd

def convert_data_to_long_df_div(species_dfs, depth):
    """Convert divergence data of RECollector to a long-form DataFrame.
    
    Parameters
    ----------
    species_dfs : dict
        Dictionary containing the dataframes for each species.
        Key: species; value: dataframe.

    depth : str
        Column of the dataframe that will be selected.

    Returns
    -------
    long_df_div : `pandas.DataFrame`
        Dataframe containing the divergence data for every species
        analyzed by RECollector. Composed of three columns: species,
        category (selected by depth), and the percentage of
        divergence for the repeat.
    """
    long_df_div = pd.DataFrame()
    for species in species_dfs:
        sp_df = species_dfs[species]
        #Create new column for species
        sp_df.insert(0, "species", species)
        selected_columns = ["species", depth, "per div"]
        long_df = sp_df.loc[:, selected_columns]
        long_df_div = pd.concat([long_df_div, long_df], axis= 0)
    long_df_div = long_df_div.reset_index(drop=True)
    return long_df_div

def get_large_dfs(file, exclude=False, transpose=False):
    """Creates DataFrame from large files.
    
    Used for RECollector's outputs.

    Parameters
    ----------
    file : file
        File containing the data (from RECollector).
        For TE count matrix also use transpose.

    exclude : bool, default: False
        If True, unknown data will be excluded.

    transpose: bool, default: False
        If True, the data will be transposed.

    Returns
    -------
    df_concat : `pandas.DataFrame`
        DataFrame with the applied changes.
    """
    df_chunk = pd.read_csv(file, header=0, index_col=0, chunksize=1000000)
    chunk_list = []
    for chunk in df_chunk:
        #For RECollector TE count matrix
        if exclude and transpose:
                chunk = chunk.T
                excluded = ["Unknown", "none", "{'none':'none'}"]
                chunk = chunk.drop(columns=excluded, errors="ignore")

        #For RECollector divergence data
        elif exclude and not transpose:
            cat_name = chunk.columns[1]
            excluded = ["Unknown", "none", "{'none':'none'}"]
            chunk = chunk[chunk[cat_name].apply(lambda x: x not in excluded)]

        #For RECollector TE count matrix
        elif not exclude and transpose:
            chunk = chunk.T

        chunk_list.append(chunk)

    df_concat = pd.concat(chunk_list)

    return df_concat

def read_chroms_file(chr_file):
    """Reads the file for chromosome filtering.
    
    Used for RECollector's chromosome filtering step.

    Parameters
    ----------
    chr_file : file
        File containing data to filter by chromosomes.
        Each row corresponds to a given species: the name
        of the species must first provided followed by
        its specified chromosomes, separated by a tab.
        For several chromosomes, they must separated by
        commas.
    
    Returns
    -------
    chrs_to_filter : dict
        Keys correspond to the names of the species, whereas
        values consist of a list of the selected chromosomes.
    """
    chrs_to_filter = {}
    for species_chrom in chr_file:
        species_chrom = species_chrom.strip("\n").split("\t")
        sp_name = species_chrom[0]
        chrs = species_chrom[1].split(",")
        chrs_to_filter[sp_name] = chrs
    return chrs_to_filter

def read_doms_file(doms_file):
    """Reads the file for domain filtering
    
    Used for RECollector's domain filtering step.

    Parameters
    ----------
    doms_file : file
        File containing two rows. The first row contains
        three tab-separated columns: domains, clades, and
        features. The second row contains the values for each
        column, separated by tabs. When creating the file, if
        one of the columns does not contain a value, use a tab.
        Values in the same column must be separated by commas.
        For the domain column, each feature must be joined by a
        colon (:), e.g. Dom1:Clade1,Dom2:Clade2.

    Returns
    -------
    domains : list
        List of domains (empty if none are provided).

    clades : list
        List of clades (empty if none are provided).

    features_dict : list of dictionaries
        List of dictionaries in which each dictionary is
        a feature (e.g. {Dom1: Clade1}), empty if none 
        are provided. 
    """
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

def read_names_file(names_file):
    """Reads the file for association of the values of its rows.
    
    Parameters
    ----------
    names_file : file
        File containing values for their association.
        Each row consists of the name of the key, and its
        value (for example, the name of the row's species 
        or the name of the group it belongs), separated by a tab.

    Returns
    -------
    row_association: dict
        Dictionary containing each row association.
    """
    row_association = {}
    for pair in names_file:
        pair = pair.strip("\n").split("\t")
        key_name = pair[0]
        value = pair[1]
        row_association[key_name] = value
    return row_association