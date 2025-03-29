import csv
import pandas as pd
import numpy as np
def parse_rescue(path):
    rows = []
    summary = []

    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if line[0][0] == "#":
                parts = line[0].split()
                dom = parts[1][:-1]  # Removes the trailing colon
                direct = int(parts[3])  # After 'DirectHits:'
                rescue = int(parts[6])  # After 'Rescued Proteins:'
                found = int(parts[12])  # Number before '('
                total = int(parts[-1][:-1])  # Last number, removing trailing ")"
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
                evalue = float(parts[1].split(":")[1])
                rounds = 1 if "1" in parts[2] else 2 if "2" in parts[2] else 0
                rows.append([sys_id, sys_len, dom_id, start, end, evalue, rounds])

    # Rescue Dataframe
    df_rows = pd.DataFrame(rows, columns=["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
    
    # Domain Summary Dataframe
    df_summary = pd.DataFrame(summary, columns=["domain", "direct hits", "rescued", "found", "total"])
    df_summary["perc found"] = df_summary["found"] / df_summary["total"]
    df_summary = df_summary.sort_values(["direct hits", "perc found"], ascending=[False, False])
    domains_priority = list(df_summary["domain"])

    filtered_domains = []  # Store the final non-overlapping domains
    
    for i, dom1 in enumerate(domains_priority):
        keep = True  # Assume we keep the domain unless an overlap is found
        
        for dom2 in filtered_domains:  # Check against already accepted domains
            if is_overlap(dom1, dom2, df_rows):
                keep = False  # Overlapping with a higher-priority domain
                break  # No need to check further
        if keep:
            filtered_domains.append(dom1)  # Keep the domain if no overlap found
    
    # Filter out domains with no overlap
    df_rows = df_rows[df_rows["subject accs."].isin(filtered_domains)]
    return df_rows

def is_overlap(dom1, dom2, df_doms):
    # Get all systems
    accessions = df_doms["query acc."].unique()
    # Check each system
    for acc in accessions:
        # Get all domains for current system
        curr = df_doms[df_doms["query acc."] == acc]
        dom1_locs = curr[curr["subject accs."] == dom1][["q. start", "q. end"]]
        dom2_locs = curr[curr["subject accs."] == dom2][["q. start", "q. end"]]
        # If any significant overlap is found
        if overlap_helper(dom1_locs, dom2_locs):
            return True
    return False

def overlap_helper(df1, df2, mutual_overlap=0.2):
    """
    Checks if any intervals in df1 overlap with any intervals in df2.

    Parameters:
    df1, df2: DataFrames with "q. start" and "q. end" columns representing intervals.

    Returns:
    True if there is at least one significant overlapping interval (mutual overlapp over longer domain >= 0.2), otherwise False.
    """
    for start1, end1 in df1.itertuples(index=False):
        for start2, end2 in df2.itertuples(index=False):
            if start1 <= end2 and start2 <= end1:  # Overlap condition
                # Overlap region divided by longer domain
                mutual = (min(end1, end2) - max(start1, start2)) / (max(end1 - start1, end2 - start2))
                if mutual >= mutual_overlap:
                    return True
    return False # no overlap


print(parse_rescue("1.A.12_rescuedDomains.tsv"))