import pandas as pd

def merge_inputs(rm_input, te_input):
    """Merges inputs from both RepeatMasker and TESorter.
    
    Repeats of RepeatMasker that are not in TESorter
    are filled with columns of the latter.

    Parameters
    ----------
    rm_input : list of dictionaries
        List of dictionaries from `read_repeatmasker_out`.
    
    te_input : list of dictionaries 
        List of dictionaries from `read_tesorter_cls_tsv`.

    Returns
    -------
    merged_inputs : `pandas.DataFrame`
    """
    #Create and concatenate dataframes
    rm_df = pd.DataFrame(rm_input)
    te_df = pd.DataFrame(te_input)
    concat_df = pd.concat([rm_df, te_df])

    #Fill empty values, domains contain a list
    tr = concat_df["domains"].isna()
    concat_df.loc[tr, "domains"] = concat_df.loc[tr, "domains"].apply(lambda x: [{"none":"none"}])
    first_cols = list(rm_df[rm_df.columns[~rm_df.columns.isin(['start','end', 'repeat'])]])
    last_cols = list(te_df[te_df.columns[~te_df.columns.isin(['start','end', 'repeat'])]])
    none_cols = [col for col in last_cols if col != "domains"]
    concat_df[none_cols] = concat_df[none_cols].fillna("none")

    #Merge rows of RM and TES belonging to the same repeat
    cols = first_cols + last_cols
    cols_order = {col:("first" if col in first_cols else "last") for col in cols} 
    merged_df = concat_df.groupby(["start", "end", "repeat"], as_index=False).aggregate(cols_order).reindex(columns=concat_df.columns)
    merged_df["id"] = merged_df["id"].astype("int32")
    merged_df = merged_df.sort_values(by=["id"]).reset_index(drop=True)

    #Convert values of numerical columns into floats/integers
    convert_dict = {
        "sw": "int32", "per div": "float32", "per del": "float32",
        "per ins": "float32","start": "int64", "end": "int64",
        "q left": "int64", "r start": "int32", "r end": "int32",
        "r left": "int32", "id": "int32", "length": "int32"
        }
    merged_df = merged_df.astype(convert_dict)

    return merged_df

def read_repeatmasker_out(input_fhand):
    """Reads the .out file from RepeatMasker.
    
    Each repeat is added as a dictionary, containing the data
    from each column as key-value pairs, to a common list for
    the file. Additionally, it separates the class/family
    column into two separate columns, and reallocates the columns
    of the start and (left) of "position in repeat" of repeats
    whose match is in the complementary strand.

    Parameters
    ----------
    input_fhand : file
        .out file from RepeatMasker.

    Returns
    -------
    rm_input : list of dictionaries
    """
    fieldnames = [
        "sw", "per div", "per del", "per ins", "seqid", "start",
        "end", "q left", "match", "repeat", "class/family", "r start",
        "r end", "r left", "id"
        ]
    read_repeats = []
    
    for _ in range(3):
        next(input_fhand)

    for line in input_fhand:
        repeat_data = {}
        line = line.strip(" ").split()

        for i in range(15):
            if i == 10:
                class_family = line[i].split("/")
                if len(class_family) == 1:
                    class_family.append("Unknown")
                cls = class_family[0]
                superfamily = class_family[1]
                repeat_data["class"] = cls
                repeat_data["superfamily"] = superfamily
            else:
                repeat_data[fieldnames[i]] = line[i]

        if repeat_data["match"] == "C":
            true_left = repeat_data["r start"]
            true_start = repeat_data["r left"]
            repeat_data["r start"] = true_start
            repeat_data["r left"] = true_left

        repeat_data["q left"] = repeat_data["q left"].strip("()")
        repeat_data["r left"] = repeat_data["r left"].strip("()")
        repeat_data["length"] = int(repeat_data["end"]) - int(repeat_data["start"])

        read_repeats.append(repeat_data)

    return read_repeats

def read_tesorter_cls_tsv(input_fhand):
    """Reads the .cls.tsv file from TESorter.
    
    Sequences that were analyzed by TESorter must be 
    previously analyzed by RepeatMasker (and have to be
    prepared with RepeatModeler). It creates a list
    of dictionaries, each dictionary corresponding to a repeat.
    It also reads data contained in the #TE column, which
    comes from RepeatMasker.

    Parameters
    ----------
    input_fhand : file
        .cls.tsv file from TESorter.

    Returns
    -------
    te_input : list of dictionaries
    """
    fieldnames = [
        "seqid", "start", "end", "repeat", "class",
        "superfamily", "tes order", "tes superfamily",
        "clade", "complete", "strand", "domains", "length"
        ]
    read_repeats = []
    next(input_fhand)

    for line in input_fhand:
        line = line.strip("\n").split("\t")

        #Classify data from the #TE column
        seqid_and_info = line[0].split(":", 1)
        seqid = seqid_and_info[0]
        info = seqid_and_info[1].split("_", 1)
        pos = info[0].split("..")
        start = pos[0]
        end = pos[1]
        length = int(end) - int(start)
        rep_info = info[1].split("#")
        rep = rep_info[0]

        class_family = rep_info[1].split("/")
        if len(class_family) == 1:
            class_family.append("Unknown")
        cls = class_family[0]
        family = class_family[1]

        #Get data from the other columns of TESorter
        tes_order = line[1]
        tes_superfamily = line[2]
        clade = line[3]
        complete = line[4]
        strand = line[5]

        #Data in the domains column will consist
        #of a list of dictionaries
        old_domains = line[6].split()
        new_domains = []
        for domain in old_domains:
            domain_clade = domain.split("|")
            if len(domain_clade) == 1:
                domain_clade.append("none")
            new_domains.append({domain_clade[0]: domain_clade[1]})

        #Create the dictionary for the repeat
        rep_data = [seqid, start, end, rep, cls,
                    family, tes_order, tes_superfamily,
                    clade, complete, strand, new_domains, length]
        repeat = dict(zip(fieldnames, rep_data))

        read_repeats.append(repeat)

    return read_repeats
