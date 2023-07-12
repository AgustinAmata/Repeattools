import numpy as np
import pandas as pd

def create_df_from_parsed_input(parsed_input):
    repeats_df = pd.DataFrame(parsed_input)
    convert_dict = {
        "sw": "int64", "per div": "float64", "per del": "float64",
        "per ins": "float64","start": "int64", "end": "int64", "r start": "int64",
        "r end": "int64", "id": "int64"
        }
    repeats_df = repeats_df.astype(convert_dict)
    return repeats_df

def filter_df_by_domain(df_to_filter, filter_by_domain, clades, special_features):
    filtered_df = df_to_filter[df_to_filter.domains.apply(lambda x: x != [{"none": "none"}])]

#Example: ["RH", "TPase"]
    if filter_by_domain:
        filtered_df = filtered_df.loc[filtered_df.domains.isin(filter_by_domain)]

#Example: ["Ale", "Ivana"]
    if clades:
        filtered_df = filtered_df.loc[filtered_df.clade.isin(clades)]

#Example: [{"RT": "Ale"}, {"RH": "Ogre"}]; function will return 
#all rows in which at least one of the given elements is found in the column
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
            print("Please select between higher_than, lower_than and equal")

    else:
        print("Please select between div, del and ins")