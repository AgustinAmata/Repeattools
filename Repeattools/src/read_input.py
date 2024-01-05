import pandas as pd
from pathlib import Path

#Get file for TE classification and process it
te_file = Path(__file__).parent.parent.absolute().joinpath("static/TEClasses_RECollector.tsv")
tes_rm_equiv_file = Path(__file__).parent.parent.absolute().joinpath("static/TES-RM_equiv.tsv")

#Data for new classification system
classifier = {}
with open(te_file) as classif_file:
    for line in classif_file:
        line = line.strip("\n")
        element_classification = line.split("\t")
        element = element_classification[0]
        classification = element_classification[1]
        classifier[element] = classification.split(";")

#Data to detect equal RepeatMasker classification values in TESorter
tes_rm_dict = {}
with open(tes_rm_equiv_file) as equiv_file:
    for line in equiv_file:
        line = line.strip("\n")
        tes_rm = line.split("\t")
        tes_cf = tes_rm[0]
        rm_equiv = tes_rm[1]
        tes_rm_dict[tes_cf] = rm_equiv

def merge_inputs(target_df, te_df, override=False):
    """Merges inputs from both RepeatMasker and TESorter.
    
    Repeats of RepeatMasker that are not in TESorter
    are filled with columns of the latter.

    It also creates a new classification system based on
    the DFAM transposable element classification to unify
    data from both RepeatMasker and TESorter. The classification
    is composed of four columns: class, subclass, superfamily
    and element. Furthermore, clade data of Copia and Gypsy
    superfamilies (from TESorter) is included during the
    processing. Unknown repeats from RepeatMasker can be
    overriden with data from TESorter.

    Parameters
    ----------
    target_df : `pandas.DataFrame`
        `pandas.DataFrame` from `read_repeatmasker_out`.
    
    te_df : `pandas.DataFrame`
        `pandas.DataFrame` from `read_tesorter_cls_tsv`.

    override : bool, default:False
        If True, unknown repeats from RepeatMasker will be
        rewritten with classification data from TESorter.
        It will also rewrite repeats which have only been
        identified as Class II transposons, given that their
        TESorter classification also belongs to this Class.

    Returns
    -------
    merged_inputs : `pandas.DataFrame`
    """
    #Select columns to merge and datatypes for the merged df
    merge_cols = ["seqid", "start", "end", "class/family", "length", "repeat"]
    cats_dict = {"class/family": "category",
                 "repeat": "category", "tes_order":"category",
                 "tes_superfamily": "category", "clade": "category"}
    target_df = target_df.merge(te_df, how="left", on=merge_cols).astype(cats_dict)
    
    #Remove repeat, start, and end columns as they are no longer necessary
    target_df.drop(["repeat", "start", "end"],inplace=True,axis=1)
    
    #Add "Unknown" to repeats of RepeatMasker that do not have
    #a match in TESorter
    none_cols = ["tes_order", "tes_superfamily", "clade"]
    for col in none_cols:
        if "Unknown" not in list(target_df[col].cat.categories):
            target_df[col] = target_df[col].cat.add_categories("Unknown").fillna("Unknown")
        else:    
            target_df[col] = target_df[col].fillna("Unknown")

    #Same as before, but with the "domains" column
    tr = target_df["domains"].isna()
    target_df.loc[tr, "domains"] = target_df.loc[tr, "domains"].apply(lambda x: [{"none":"none"}])

    #Create new classification for the data
    classif_cols = ["class", "subclass", "superfamily", "element"]
    cl_fam_col = target_df["class/family"].astype("str")
    for i, col in enumerate(classif_cols):
        if i==3: #Add "clade" data from Copia and Gypsy to the classification
            target_df[col] = cl_fam_col.apply(lambda x: classifier[x][i])
            #Get Copia and Gypsy rows that match with their TESorter data and add the "clade" data
            copia_gypsy_values = (target_df["superfamily"].astype(str).str.contains("Copia|Gypsy")) & (target_df["element"].astype(str)=="Unknown") & (target_df["tes_superfamily"].astype(str).str.contains("Copia|Gypsy"))
            target_df.loc[copia_gypsy_values, col] = target_df.loc[copia_gypsy_values, "clade"].astype("str")
            target_df[col] = target_df[col].astype("category")
        else:
            target_df[col] = cl_fam_col.apply(lambda x: classifier[x][i]).astype("category")

    #Override Unknown and Class_II unknown (only if TESorter data matches) RepeatMasker rows
    if override:
        target_df["tes_classif"] = pd.Series(target_df["tes_order"].astype("str") + "/" + target_df["tes_superfamily"].astype("str")).astype("category")
        target_df["tes_classif"] = target_df["tes_classif"].cat.rename_categories(tes_rm_dict)
        class_unknown_values = target_df["class"]=="Unknown"
        dna_unknown_values = (target_df["class"].astype(str).str.contains("Class_II")) & (target_df["subclass"].astype(str)=="Unknown") & (target_df["tes_order"].astype(str).str.contains("TIR|Helitron|Maverick"))
        class_and_dna_values = class_unknown_values + dna_unknown_values
        for i, col in enumerate(classif_cols):
            target_df[col] = target_df[col].astype("str")
            target_df.loc[class_and_dna_values, col] = target_df.loc[class_and_dna_values, "tes_classif"].astype("str").apply(lambda x: classifier[x][i])
            if i == 3: #Add "clade" data from Copia and Gypsy to the classification
                #Get Copia and Gypsy rows that match with their TESorter data and add the "clade" data
                copia_gypsy_values = (target_df["superfamily"].astype(str).str.contains("Copia|Gypsy")) & (target_df["element"].astype(str)=="Unknown") & (target_df["tes_superfamily"].astype(str).str.contains("Copia|Gypsy"))
                target_df.loc[copia_gypsy_values, col] = target_df.loc[copia_gypsy_values, "clade"].astype("str")
            target_df[col] = target_df[col].astype("category")

    #Finally, remove "seqid", "tes_classif", and "class/family" columns
    target_df.drop(["seqid", "tes_classif", "class/family"],inplace=True,axis=1)

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
    # def create_cf_list(line):
    #     class_family = line.split("/")
    #     if len(class_family) == 1:
    #         class_family.append("Unknown")
    #     return class_family

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
        "end": "int32", "repeat": "category", "class/family": "category"
        }

    #Create the DataFrame, and
    #add a repeat length column
    rm_input = pd.read_csv(input_fhand, skiprows=2, header=None,
                           names=fieldnames, sep="\s+", comment="*",
                           usecols=usable_cols, dtype=convert_dict)
    rm_input["length"] = rm_input["end"] - rm_input["start"]

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

    #Initial columns of the the file, columns that will be
    #used and their datatypes    
    fieldnames = [
        "#TE", "tes_order", "tes_superfamily", "clade",
        "complete", "strand", "domains"
        ]
    usable_cols = [
        "#TE", "tes_order", "tes_superfamily", "clade",
        "domains"
        ]
    convert_dict = {"tes_order": "category",
                    "tes_superfamily": "category",
                    "clade": "category"}
    
    #Create the DataFrame
    te_input = pd.read_csv(input_fhand, delimiter="\t", header=0,
                           names=fieldnames, usecols=usable_cols,
                           dtype=convert_dict)

    #Separate data of #TE column
    new_cols = ["seqid", "start", "end", "repeat", "class/family"]
    te_input[new_cols] = te_input.pop("#TE").str.extract("(.*):(\d*)..(\d*)_{1}(.*)#(.*)", expand=True)
    te_input["seqid"] = te_input["seqid"].astype("str")
    te_input["repeat"] = te_input["repeat"].astype("category")
    te_input["class/family"] = te_input["class/family"].astype("category")
    te_input[["start","end"]] = te_input[["start","end"]].astype("int32")

    #Capitalize "unknown" data of "tes_order", "tes_superfamily" and "clade"
    for col in ["tes_order", "tes_superfamily", "clade"]:
        te_input[col] = te_input[col].cat.rename_categories({"unknown": "Unknown"})

    #add a repeat length column and modify the domains column
    te_input["domains"] = te_input["domains"].apply(doms_dict)

    te_input["length"] = te_input["end"] - te_input["start"]

    return te_input