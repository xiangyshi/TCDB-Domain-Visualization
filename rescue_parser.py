import csv
import pandas as pd

def parse_rescue(path):
    rows = []

    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if line[0][0] == "#":
                continue
            print(line)
            sys_id = line[1].split(":")[0]
            sys_len = line[1].split(":")[1]
            for domain in line[2:]:
                parts = domain.split("|")
                dom_id = parts[0]
                pos = parts[1].split(":")[0]
                if pos == "Nohit":
                    continue
                start = pos.split("-")[0]
                end =pos.split("-")[1]
                evalue = parts[1].split(":")[1]
                rounds = 0
                if "1" in parts[2]:
                    rounds = 1
                elif "2" in parts[2]:
                    rounds = 2
                rows.append([sys_id, sys_len, dom_id, start, end, evalue, rounds])

    with open("temp.csv", "w") as out:
        writer = csv.writer(out, delimiter=",")
        writer.writerow(["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
        writer.writerows(rows)

    return pd.read_csv("temp.csv")
        
