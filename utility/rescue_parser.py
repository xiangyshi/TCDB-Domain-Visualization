import csv
import os
import pandas as pd
import numpy as np
from utility.util import is_overlap

def parse_rescue(path):
    """
    Parse a single rescue file and return its domain data.
    
    Args:
        path (str): Path to a single rescue file
        
    Returns:
        pd.DataFrame: DataFrame containing domain data from the rescue file
    """
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
    df_rows["family"] = df_rows["query acc."].apply(lambda x: '.'.join(x.split(".")[:3]))
    
    return df_rows

def clean_rescue(folder_path, target_fam_ids=None):
    """
    Process all rescue files in a folder and return a concatenated DataFrame.
    
    Args:
        folder_path (str): Path to folder containing rescue files
        target_fam_ids (list, optional): List of family IDs to process. If None, process all families.
        
    Returns:
        pd.DataFrame: DataFrame with all domain data from rescue files
    """
    if not os.path.isdir(folder_path):
        raise ValueError(f"The specified path is not a directory: {folder_path}")
    
    all_dfs = []
    
    # Process each rescue file in the folder
    for file in os.listdir(folder_path):
        if file.endswith("_rescuedDomains.tsv"):
            fam_id = file.split("_")[0]  # Extract family ID from filename
            
            # Skip if not in target families
            if target_fam_ids and fam_id not in target_fam_ids:
                continue
                
            try:
                file_path = os.path.join(folder_path, file)
                df = parse_rescue(file_path)
                if not df.empty:
                    all_dfs.append(df)
            except Exception as e:
                print(f"Error processing {file}: {e}")
                continue
    
    # Concatenate all DataFrames
    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    else:
        return pd.DataFrame()