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
    # Validate that the input path is a file
    if not os.path.isfile(path):
        raise ValueError(f"The specified path is not a file: {path}")
    
    # Initialize data structures to store parsed information
    rows = []  # Will store individual domain hits
    summary = []  # Will store domain summary statistics
    significant_domain_hits = {}  # Tracks count of significant hits per domain
    
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            # Process header lines that start with '#'
            if line[0][0] == "#":
                # Example line: '# CDD166458:  DirectHits:  9    Rescued Proteins:  2    Prots with Domain in 1.A.12:  11 (100.0% from a total of 11)'
                dom = line[0].split(":")[0][2:]  # Extract domain ID (removing '# ' prefix)
                found = 0  # Initialize found count
                total = int(line[0].split("of")[1][:-1].strip())  # Extract total protein count
                summary.append([dom, found, total])
                continue
            
            # Process data lines
            sys_id, sys_len = line[1].split(":")  # Extract system ID and length
            sys_len = int(sys_len)
            
            # Process each domain hit in the line
            for domain in line[2:]:
                print()
                print(domain)
                parts = domain.split("|")  # Split domain information
                dom_id = parts[0]  # Domain ID
                
                # Skip if no hit was found
                if len(parts) < 3:
                    continue
                
                round = parts[-1]

                for i in range(1, len(parts) - 1):
                    # Extract start and end positions
                    pos = parts[i].split(":")[0]  # Position information
                    start, end = pos.split("-")
                    start = int(start)
                    end = int(end)

                    # Track significant hits (DirectHit or Rescued1)
                    
                    print(dom_id, pos, start, end)
                    if round in ["DirectHit", "Rescued1"]:
                        if dom_id in significant_domain_hits:
                            significant_domain_hits[dom_id] += 1
                        else:
                            significant_domain_hits[dom_id] = 1

                    # Extract e-value and determine rescue round
                    evalue = float(parts[1].split(":")[1])
                    rounds = 1 if "1" in parts[2] else 2 if "2" in parts[2] else 0
                    rows.append([sys_id, sys_len, dom_id, start, end, evalue, rounds])
    
    # Update summary with actual found counts
    for row in summary:
        try:
            row[1] = significant_domain_hits[row[0]]
        except KeyError:
            continue

    # Create DataFrame from parsed domain hits
    df_rows = pd.DataFrame(rows, columns=["query acc.", "query length", "subject accs.", "q. start", "q. end", "evalue", "rescue round"])
    
    # Create and process summary DataFrame
    df_summary = pd.DataFrame(summary, columns=["domain", "found", "total"])
    df_summary["perc found"] = df_summary["found"] / df_summary["total"]  # Calculate percentage of found domains
    # Filter domains with at least 80% found rate and sort by percentage
    df_summary = df_summary[df_summary["perc found"] >= 0.8].sort_values("perc found", ascending=False)

    # Get list of filtered domains
    filtered_domains = list(df_summary["domain"])
    # Filter rows to only include domains that meet the threshold
    df_rows = df_rows[df_rows["subject accs."].isin(filtered_domains)]
    # Add family information by extracting first three parts of the query accession
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
    # Validate that the input path is a directory
    if not os.path.isdir(folder_path):
        raise ValueError(f"The specified path is not a directory: {folder_path}")
    
    all_dfs = []  # List to store DataFrames from each processed file
    print(f"Processing rescue files in {folder_path}, target_fam_ids: {target_fam_ids}")
    
    # Process each rescue file in the folder
    for file in os.listdir(folder_path):
        if file.endswith("_rescuedDomains.tsv"):
            fam_id = file.split("_")[0]  # Extract family ID from filename
            # Skip if not in target families (if specified)
            if target_fam_ids and (fam_id not in target_fam_ids):
                continue
            print(f"Processing {file}")
            try:
                file_path = os.path.join(folder_path, file)
                print(f"Processing {file_path}")
                df = parse_rescue(file_path)
                if not df.empty:
                    all_dfs.append(df)
            except Exception as e:
                print(f"Error processing {file}: {e}")
                continue
    
    # Concatenate all DataFrames if any were successfully processed
    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    else:
        return pd.DataFrame()  # Return empty DataFrame if no files were processed
    

if __name__ == "__main__":
    # Example usage: parse a single rescue file
    print(parse_rescue("./rescued/2.A.123_rescuedDomains.tsv"))