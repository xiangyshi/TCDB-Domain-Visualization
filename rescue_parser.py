import csv
import os
import pandas as pd
import numpy as np
from util import is_overlap

def parse_rescue(path):
    # Validate path
    if not os.path.exists(path):
        parse_path = f"rescued/{path}_rescuedDomains.tsv"
        if not os.path.exists(parse_path):
            raise FileNotFoundError(f"Rescue file not found at path: {path} or {parse_path}")
        path = parse_path
    
    if not os.path.isfile(path):
        raise ValueError(f"The specified path is not a file: {path}")
    
    rows = []
    summary = []
    significant_domain_hits = {}
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if line[0][0] == "#":
                # Line = ['# CDD166458:  DirectHits:  9    Rescued Proteins:  2    Prots with Domain in 1.A.12:  11 (100.0% from a total of 11)']
                dom = line[0].split(":")[0][2:] # domain before first colon
                found = 0
                total = int(line[0].split("of")[1][:-1].strip()) # total proteins after "of"
                summary.append([dom, found, total])
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
        try:
            row[1] = significant_domain_hits[row[0]]
        except KeyError:
            continue

    # Rescue Dataframe
    df_rows = pd.DataFrame(rows, columns=["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
    
    # Domain Summary Dataframe
    df_summary = pd.DataFrame(summary, columns=["domain", "found", "total"])
    df_summary["perc found"] = df_summary["found"] / df_summary["total"]
    df_summary = df_summary[df_summary["perc found"] >= 0.8].sort_values("perc found", ascending=False)

    filtered_domains = list(df_summary["domain"])
    # Filter out domains with no overlap
    df_rows = df_rows[df_rows["subject accs."].isin(filtered_domains)]
    return df_rows