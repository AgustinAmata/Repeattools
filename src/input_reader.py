from csv import DictReader

def read_ltr_retriever_list(input_fhand):
    read_repeats = []
    for line in DictReader(input_fhand, delimiter= " "):
        read_repeats.append(line)

    return read_repeats

#for gff3 files of LTR_Retriever v2.8.7+
def read_ltr_retriever_gff3(input_fhand):
    read_repeats = []
    for line in DictReader(
            input_fhand, fieldnames= ["seqid", "source",
            "repeat_class/superfamily", "start", "end",
            "sw_score", "strand", "phase", "attributes"],
            delimiter = "\t"
            ):

        if line["seqid"].startswith("##"):
            continue
        
        read_repeats.append(line)

    return read_repeats

def read_repeatmasker(input_fhand):
    fieldnames = [
        "sw", "per div", "per del", "per ins", "query", "q begin",
        "q end", "q left", "repeat", "class/family", "r begin",
        "r end", "r left", "id"
        ]
    read_repeats = []
    for line in input_fhand:
        repeat_data = {}
        line = line.strip(" ").split()

        for i in range(len(line)):
            repeat_data[fieldnames[i]] = line[i]

        read_repeats.append(repeat_data)

    return read_repeats