import pandas as pd

def count_tes(input_df, species_name, col="superfamily"):
    """Counts each element of the selected column.

    Parameters
    ----------
    input_df : `pandas.DataFrame`
        Dataframe in which rows correspond to repeats and some
        of the columns consist of the available categories.

    species_name : str
        Name of the species to name the resulting Series.

    col : str, default: 'superfamily'
        Column that will be selected and counted.
        
    Returns
    -------
    counted_tes : `pandas.Series`
        Indexes are named after counted elements and the Series
        is named after the species it came from, so that it
        can be used as the name of the column when combining
        different Series.
    """

    if col == "domains":
        input_df["domains"] = input_df['domains'].apply(lambda x: ','.join(map(str, x)))

    counted_tes = input_df.value_counts(col).rename(species_name).astype("int32")
    return counted_tes

def create_te_count_matrix(list_of_inputs):
    """Combines the series from count_tes() into a dataframe.
    

    Parameters
    ----------
    list_of_inputs : list
        List containing all the `pandas.Series` to concatenate.

    Returns
    -------
    te_count_matrix : `pandas.DataFrame`
        Columns: name of the species; indexes: element. In case
        some element was not present in a species, it is filled
        with a 0. All numbers are converted into integers,
        as some Series can be added as floats.
    """
    te_count_matrix = pd.DataFrame()
    for input in list_of_inputs:
        te_count_matrix = pd.concat([te_count_matrix, input], axis=1)
    te_count_matrix = te_count_matrix.fillna(0).astype("int32")
    return te_count_matrix

def filter_df_by_domain(df_to_filter, doms, clades, special_features):
    """Filters a dataframe to remove repeats that do not contain data on their domains.
    
    It can also filter by specific domains, clades and domain:clade features.

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    doms : list, optional
        List of domains to filter (include).

    clades : list, optional
        List of domains to filter (include).

    special_features : list of dictionaries, optional
        List of dictionaries in which each dictionary specifies
        a different domain:clade feature.

    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.
    """
    #Filtering is based on True/False series,
    #checking if there is at least one of the given elements in any row
    filtered_df = df_to_filter[df_to_filter.domains.apply(lambda x: x != [{"none": "none"}])]

    if doms:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.get(dom) for dom in doms) for feat in x))]

    if clades:
        filtered_df = filtered_df[filtered_df.clade.isin(clades)]

    if special_features:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.items() == dom.items() for feat in special_features) for dom in x))]

    return filtered_df

def filter_df_by_length(df_to_filter, length):
    """Filters the dataframe by eliminating rows whose length is lower than the given.

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    length : int
        Minimum length to remain in the dataframe.

    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.    
    """
    filtered_df = df_to_filter[df_to_filter["length"] >= length]

    return filtered_df

def filter_df_by_percentages(df_to_filter, percentage="div=20.0", mode="lower_than"):
    """Filters the dataframe by percentages and mode

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    percentage : str, default: 'div=20.0', optional
        Percentage to filter the dataframe. The number must
        be preceded by either 'div', 'del' or 'ins', and joined
        by a an equal sign (=).
    
    mode : str, default: 'lower_than', optional
        Mode to filter data to only include repeats 
        lower, higher or equal to the threshold
    
    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.    
    """
    name_and_perc = percentage.split("=")
    name = "per " + name_and_perc[0]
    perc = float(name_and_perc[1])
    names = ["per div", "per del", "per ins"]      
    modes = ["higher_than", "lower_than", "equal"]

    if name in names:
        if mode in modes:
            if mode == "higher_than":
                mode_filter = df_to_filter[name] >= perc
                filtered_df = df_to_filter[mode_filter]
            elif mode == "lower_than":
                mode_filter = df_to_filter[name] <= perc
                filtered_df = df_to_filter[mode_filter]
            elif mode == "equal":
                mode_filter = df_to_filter[name] == perc
                filtered_df = df_to_filter[mode_filter]

            return filtered_df