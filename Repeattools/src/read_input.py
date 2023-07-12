from csv import DictReader

def merge_inputs(rm_input, te_input):
    merged_inputs = []
    for repeat in rm_input:
        rep_start = repeat["start"]
        rep_end = repeat["end"]
        rep_name = repeat["repeat"]
        matching_rep = next((rep for rep in te_input if rep["start"] == rep_start and rep["end"] == rep_end and rep["repeat"] == rep_name), None)

        if not matching_rep:
            repeat["domains"] = [{"none": "none"}]
            key_names = ["tes order", "tes superfamily", "complete", "strand", "clade"]
            for key in key_names:
                repeat[key] = "none"
            merged_rep = repeat

        else:
            merged_rep = repeat | matching_rep

        merged_inputs.append(merged_rep)

    return merged_inputs

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

            if repeat_data["match"] == "C":
                true_left = repeat_data["r start"]
                true_start = repeat_data["r left"]
                repeat_data["r start"] = true_start
                repeat_data["r left"] = true_left

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
        new_domains = []

        for domain in old_domains:
            domain_clade = domain.split("|")
            if len(domain_clade) == 1:
                domain_clade.append("none")
            new_domains.append({domain_clade[0]: domain_clade[1]})
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