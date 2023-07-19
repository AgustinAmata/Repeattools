import pandas as pd

def create_df_from_parsed_input(parsed_input):
    repeats_df = pd.DataFrame(parsed_input)
    convert_dict = {
        "sw": "int64", "per div": "float64", "per del": "float64",
        "per ins": "float64","start": "int64", "end": "int64",
        "r start": "int64", "r end": "int64", "id": "int64"
        }
    repeats_df = repeats_df.astype(convert_dict)
    return repeats_df

def count_tes(input_df, species_name, depth="superfamily"):
    counted_tes = input_df.value_counts(depth).rename(species_name)
    return counted_tes

def create_te_count_matrix(list_of_inputs):
    te_count_matrix = pd.DataFrame()
    for input in list_of_inputs:
        te_count_matrix = pd.concat([te_count_matrix, input], axis=1)
    te_count_matrix = te_count_matrix.fillna(0)
    return te_count_matrix.astype("int64")

def filter_df_by_chromosomes(df_to_filter, chromosomes, exclude=False):
    if exclude:
        filtered_df = df_to_filter[~df_to_filter.seqid.apply(lambda x: x in chromosomes)]

    else:
        filtered_df = df_to_filter[df_to_filter.seqid.apply(lambda x: x in chromosomes)]

    return filtered_df

def filter_df_by_domain(df_to_filter, doms, clades, special_features):
    filtered_df = df_to_filter[df_to_filter.domains.apply(lambda x: x != [{"none": "none"}])]

    if doms:
        filtered_df = filtered_df.loc[filtered_df.domains.apply(lambda x: any(any(feat.get(dom) for dom in doms) for feat in x))]

    if clades:
        filtered_df = filtered_df.loc[filtered_df.clade.isin(clades)]

    if special_features:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.items() == dom.items() for feat in special_features) for dom in x))]

    return filtered_df

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