import csv
import pandas as pd

def parse_rescue(path):
    rows = []
    summary = []

    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if line[0][0] == "#":
                parts = line[0].split()
                print(parts)
                dom = parts[1][:-1]  # Removes the trailing colon
                direct = int(parts[3])  # After 'DirectHits:'
                rescue = int(parts[6])  # After 'Rescued Proteins:'
                found = int(parts[12])  # Number before '('
                total = int(parts[-1][:-1])  # Last number, removing trailing ")"
                summary.append([dom, direct, rescue, found, total])
                continue

            sys_id, sys_len = line[1].split(":")
            for domain in line[2:]:
                parts = domain.split("|")
                dom_id = parts[0]
                pos = parts[1].split(":")[0]
                if pos == "Nohit":
                    continue
                start, end = pos.split("-")
                evalue = parts[1].split(":")[1]
                rounds = 1 if "1" in parts[2] else 2 if "2" in parts[2] else 0
                rows.append([sys_id, sys_len, dom_id, start, end, evalue, rounds])

    df_rows = pd.DataFrame(rows, columns=["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
    df_summary = pd.DataFrame(summary, columns=["domain", "direct hits", "rescued", "found", "total"])
    
    return df_rows, df_summary

print(parse_rescue("1.A.12_rescuedDomains.tsv"))