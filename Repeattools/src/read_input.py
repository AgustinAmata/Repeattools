import pandas as pd

def merge_inputs(target_df, te_df):
    """Merges inputs from both RepeatMasker and TESorter.
    
    Repeats of RepeatMasker that are not in TESorter
    are filled with columns of the latter.

    Parameters
    ----------
    target_df : `pandas.DataFrame`
        `pandas.DataFrame` from `read_repeatmasker_out`.
    
    te_df : `pandas.DataFrame`
        `pandas.DataFrame` from `read_tesorter_cls_tsv`.

    Returns
    -------
    merged_inputs : `pandas.DataFrame`
    """
    #Select columns to merge and datatypes for the merged df
    merge_cols = ["seqid", "start", "end", "class", "superfamily", "length", "repeat"]
    cats_dict = {"class": "category", "superfamily": "category",
                 "repeat": "category", "tes order":"category",
                 "tes superfamily": "category", "clade": "category"}
    target_df = target_df.merge(te_df, how="left", on=merge_cols).astype(cats_dict)

    #Add "Unknown" to repeats of RepeatMasker that do not have
    #a match in TESorter
    none_cols = ["tes order", "tes superfamily", "clade"]
    for col in none_cols:
        if "Unknown" not in list(target_df[col].cat.categories):
            target_df[col] = target_df[col].cat.add_categories("Unknown").fillna("Unknown")
        else:    
            target_df[col] = target_df[col].fillna("Unknown")

    #Same as before, but with the "domains" column
    tr = target_df["domains"].isna()
    target_df.loc[tr, "domains"] = target_df.loc[tr, "domains"].apply(lambda x: [{"none":"none"}])

    return target_df

def read_repeatmasker_out(input_fhand):
    """Reads the .out file from RepeatMasker.
    
    It creates a `pandas.DataFrame`. Additionally,
    it separates the class/family column into two
    separate columns.

    Parameters
    ----------
    input_fhand : file
        .out file from RepeatMasker.

    Returns
    -------
    rm_input : `pandas.DataFrame`
    """
    #Function to prepare for separation of the class/family column
    def create_cf_list(line):
        class_family = line.split("/")
        if len(class_family) == 1:
            if class_family[0] == "Unknown":
                class_family.append("Unknown")
            else:
                class_family.append(class_family[0])
        return class_family

    #Initial columns of the the file, columns that will be
    #used and their datatypes
    fieldnames = [
        "sw", "per div", "per del", "per ins", "seqid", "start",
        "end", "q left", "match", "repeat", "class/family", "r start",
        "r end", "r left", "id"
        ]
    usable_cols = [
        "per div", "per del", "per ins", "seqid", "start",
        "end", "repeat", "class/family"
        ]
    convert_dict = {
        "per div": "float16", "per del": "float16",
        "per ins": "float16", "seqid": "str", "start": "int32",
        "end": "int32", "repeat": "category"
        }

    #Create the DataFrame, separate the class/family column and
    #add a repeat length column
    rm_input = pd.read_csv(input_fhand, skiprows=2, header=None,
                           names=fieldnames, sep="\s+", comment="*",
                           usecols=usable_cols, dtype=convert_dict)
    rm_input["length"] = rm_input["end"] - rm_input["start"]
    rm_input["class/family"] = rm_input["class/family"].apply(create_cf_list)
    rm_input[["class", "superfamily"]] = pd.DataFrame(rm_input.pop("class/family").tolist(),
                                                      index=rm_input.index, dtype="category")

    return rm_input

def read_tesorter_cls_tsv(input_fhand):
    """Reads the .cls.tsv file from TESorter.
    
    Sequences that were analyzed by TESorter must be 
    previously analyzed by RepeatMasker (and have to be
    prepared with RepeatModeler). It creates a `pandas.DataFrame`.
    It also reads data contained in the #TE column, which
    comes from RepeatMasker.

    Parameters
    ----------
    input_fhand : file
        .cls.tsv file from TESorter.

    Returns
    -------
    te_input : `pandas.DataFrame`
    """
    #Function to create list of dictionaries for each repeat domains
    def doms_dict(dom):
        old_dom = dom.split()
        new_doms = []
        for dom in old_dom:
            domain_clade = dom.split("|")
            if len(domain_clade) == 1:
                domain_clade.append("none")
            new_doms.append({domain_clade[0]:domain_clade[1]})
        return new_doms

    #Function to prepare for separation of the class/family column
    def create_cf_list(line):
        class_family = line.split("/")
        if len(class_family) == 1:
            class_family.append("Unknown")
        return class_family
    
    #Initial columns of the the file, columns that will be
    #used and their datatypes    
    fieldnames = [
        "#TE", "tes order", "tes superfamily", "clade",
        "complete", "strand", "domains"
        ]
    usable_cols = [
        "#TE", "tes order", "tes superfamily", "clade",
        "domains"
        ]
    convert_dict = {"tes order": "category",
                    "tes superfamily": "category",
                    "clade": "category"}
    
    #Create the DataFrame
    te_input = pd.read_csv(input_fhand, delimiter="\t", header=0,
                           names=fieldnames, usecols=usable_cols,
                           dtype=convert_dict)

    #Separate data of #TE column
    new_cols = ["seqid", "start", "end", "repeat", "class/superfamily"]
    te_input[new_cols] = te_input.pop("#TE").str.extract("(.*):(\d*)..(\d*)_{1}(.*)#(.*)", expand=True)
    te_input["repeat"] = te_input["repeat"].astype("category")
    te_input[["start","end"]] = te_input[["start","end"]].astype("int32")

    #Capitalize "unknown" data of "tes order", "tes superfamily" and "clade"
    for col in ["tes order", "tes superfamily", "clade"]:
        te_input[col] = te_input[col].cat.rename_categories({"unknown": "Unknown"})

    #Separate the class/family column,
    #add a repeat length column and modify the domains column
    te_input["domains"] = te_input["domains"].apply(doms_dict)
    te_input["class/superfamily"] = te_input["class/superfamily"].apply(create_cf_list)
    te_input[["class", "superfamily"]] = pd.DataFrame(te_input.pop("class/superfamily").tolist(),
                                                      index=te_input.index, dtype="category")
    te_input["length"] = te_input["end"] - te_input["start"]

    return te_input