from csv import DictReader

def read_ltr_retriever_list(input_fhand):
    read_repeats = []
    for repeat in DictReader(input_fhand, delimiter= " "):
        repeat = {k.lower():v for (k, v) in repeat.items()}

        seqid_and_pos = repeat["#ltr_loc"].split(":")
        seqid = seqid_and_pos[0]
        pos = seqid_and_pos[1].split("..")
        start = pos[0]
        end = pos[1]

        repeat["seqid"] = seqid
        repeat["start"] = start
        repeat["end"] = end
        repeat.pop("#ltr_loc")
        read_repeats.append(repeat)

    return read_repeats

def read_repeatmasker_out(input_fhand):
    fieldnames = [
        "sw", "per div", "per del", "per ins", "seqid", "start",
        "end", "q left", "match", "repeat", "class/family", "r start",
        "r end", "r left", "id"
        ]
    read_repeats = []
    
    for line in input_fhand:
        repeat_data = {}
        line = line.strip(" ").split()

        if not line:
            continue

        elif line[0].startswith(("SW", "score")):
            continue

        else:
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

            read_repeats.append(repeat_data)

    return read_repeats

def read_tesorter_cls_tsv(input_fhand):
    read_repeats = []
    for repeat in DictReader(input_fhand, delimiter= "\t"):
        repeat = {k.lower():v for (k, v) in repeat.items()}

        seqid_and_info = repeat["#te"].split(":", 1)
        seqid = seqid_and_info[0]
        info = seqid_and_info[1].split("_", 1)
        pos = info[0].split("..")
        start = pos[0]
        end = pos[1]
        rep_info = info[1].split("#")
        rep = rep_info[0]
        class_family = rep_info[1].split("/")

        if len(class_family) == 1:
            class_family.append("Unknown")

        cls = class_family[0]
        family = class_family[1]

        old_domains = repeat["domains"].split()
        new_domains = {}

        for domain in old_domains:
            domain_clade = domain.split("|")
            if len(domain_clade) == 1:
                domain_clade.append("none")
            new_domains[domain_clade[0]] = domain_clade[1]
        repeat["domains"] = new_domains

        repeat["tes order"] = repeat["order"]
        repeat["tes superfamily"] = repeat["superfamily"]
        repeat.pop("order")

        key_names = ["seqid", "start", "end", "repeat", "class", "superfamily"]
        values = [seqid, start, end, rep, cls, family]

        for i in range(6):
            repeat[key_names[i]] = values[i]

        repeat.pop("#te")
        read_repeats.append(repeat)

    return read_repeats