import csv
import pandas as pd
import numpy as np
from util import is_overlap

def parse_rescue(path):
    rows = []
    summary = []
    significant_domain_hits = {}
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            print("processing line: ", line)
            if line[0][0] == "#":
                dom = line[0].split(":")[0][2:]  # line[0] = '# CDD238744:  DirectHits:  13'
                direct = int(line[0].split(":")[2].strip())  # line[0] = '# CDD238744:  DirectHits:  13'
                rescue = int(line[1].split(":")[1].strip())  # line[1] = 'Rescued Proteins:  3'
                found = 0
                total = int(line[2].split("of")[1][:-1].strip())  # Last number, removing trailing ")"
                summary.append([dom, direct, rescue, found, total])
                continue
            
            
            sys_id, sys_len = line[1].split(":")
            sys_len = int(sys_len)
            for domain in line[2:]:
                parts = domain.split("|")
                dom_id = parts[0]
                pos = parts[1].split(":")[0]
                if pos == "Nohit":
                    continue
                start, end = pos.split("-")
                start = int(start)
                end = int(end)

                pos = parts[2]
                if pos in ["DirectHit", "Rescued1"]:
                    if dom_id in significant_domain_hits:
                        significant_domain_hits[dom_id] += 1
                    else:
                        significant_domain_hits[dom_id] = 1

                evalue = float(parts[1].split(":")[1])
                rounds = 1 if "1" in parts[2] else 2 if "2" in parts[2] else 0
                rows.append([sys_id, sys_len, dom_id, start, end, evalue, rounds])
    for row in summary:
        row[3] = significant_domain_hits[row[0]]

    # Rescue Dataframe
    df_rows = pd.DataFrame(rows, columns=["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
    
    # Domain Summary Dataframe
    df_summary = pd.DataFrame(summary, columns=["domain", "direct hits", "rescued", "found", "total"])
    df_summary["perc found"] = df_summary["found"] / df_summary["total"]
    df_summary = df_summary[df_summary["perc found"] >= 0.8].sort_values(["direct hits", "perc found"], ascending=[False, False])

    filtered_domains = list(df_summary["domain"])
    # Filter out domains with no overlap
    df_rows = df_rows[df_rows["subject accs."].isin(filtered_domains)]
    return df_rows