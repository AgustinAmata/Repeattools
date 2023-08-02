import pandas as pd

#Creates a Pandas dataframe from the merged input
def create_df_from_parsed_input(parsed_input):
    repeats_df = pd.DataFrame(parsed_input)
    convert_dict = {
        "sw": "int32", "per div": "float32", "per del": "float32",
        "per ins": "float32","start": "int64", "end": "int64",
        "q left": "int64", "r start": "int32", "r end": "int32",
        "r left": "int32", "id": "int32", "length": "int32"
        }
    repeats_df = repeats_df.astype(convert_dict)
    return repeats_df

#Counts each element of the selected column and returns it
#as a Pandas series in which indexes are named after counted elements
#and the series is named after the species it came from, so that it
#can be used as the name of the column when combining different
#series
def count_tes(input_df, species_name, depth="superfamily"):
    if depth == "domains":
        input_df["domains"] = input_df['domains'].apply(lambda x: ','.join(map(str, x)))

    counted_tes = input_df.value_counts(depth).rename(species_name)
    return counted_tes

#Combines the series from count_tes() into a dataframe (columns:
#name of the species; indexes: element). In case some element was not
#present in a species, it fills it with a 0. Finally, it converts
#all numbers into integers, as some series can be added as floats
def create_te_count_matrix(list_of_inputs):
    te_count_matrix = pd.DataFrame()
    for input in list_of_inputs:
        te_count_matrix = pd.concat([te_count_matrix, input], axis=1)
    te_count_matrix = te_count_matrix.fillna(0)
    return te_count_matrix.astype("int32")

#Filters a dataframe based on a list of chromosomes so that it only contains
#(or excludes) those chromosomes
def filter_df_by_chromosomes(df_to_filter, chromosomes, exclude=False):
    if exclude:
        filtered_df = df_to_filter[~df_to_filter.seqid.apply(lambda x: x in chromosomes)]

    else:
        filtered_df = df_to_filter[df_to_filter.seqid.apply(lambda x: x in chromosomes)]

    return filtered_df

#Filters a dataframe to remove repeats that do not contain data on
#their domains. It can also filter by specific domains, clades and
#domain:clade features. Filtering is based on True/False series,
#checking if there is at least one of the given elements in any row
def filter_df_by_domain(df_to_filter, doms, clades, special_features):
    filtered_df = df_to_filter[df_to_filter.domains.apply(lambda x: x != [{"none": "none"}])]

    if doms:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.get(dom) for dom in doms) for feat in x))]

    if clades:
        filtered_df = filtered_df[filtered_df.clade.isin(clades)]

    if special_features:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.items() == dom.items() for feat in special_features) for dom in x))]

    return filtered_df

#Filters the dataframe by eliminating rows whose length
#is lower than the given
def filter_df_by_length(df_to_filter, length):
    filtered_df = df_to_filter[df_to_filter["length"] >= length]

    return filtered_df

#Filters the dataframe by percentages
def filter_df_by_percentages(df_to_filter, percentage="div=20.0", mode="lower_than"):
    name_and_perc = percentage.split("=")
    name = "per " + name_and_perc[0]
    perc = float(name_and_perc[1])
    names = ["per div", "per del", "per ins"]

    higher = df_to_filter[name] >= perc
    lower = df_to_filter[name] <= perc
    equal = df_to_filter[name] == perc        
    modes = {"higher_than": higher, "lower_than": lower, "equal": equal}

    if name in names:
        if mode in modes:
            filtered_df = df_to_filter[modes[mode]]
            return filtered_df

    else:
        print("Please select between div, del and ins")