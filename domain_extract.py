'''
Filename: domain_extract.py
Author: Xiangyu(Leo) Shi
Created: Aug 2024

This script processes protein domains using rpsblast and TCDB-specific proteins.
It analyzes domain data from CDD files and generates domain architecture analysis
for specified protein families.
'''

from util import *
from config import *
import argparse
import numpy as np
import pandas as pd
import csv
from Family import *

def parse_arguments():
    """
    Parses and validates command-line arguments.

    Returns:
        tuple: A tuple containing:
            - args (argparse.Namespace): Parsed command-line arguments
            - isFile (bool): Flag indicating if the TCID input is a file
    
    Raises:
        ArgumentError: If the input file doesn't exist
    """
    parser = argparse.ArgumentParser(
        description="Process protein domains using rpsblast and TCDB-specific proteins."
    )

    parser.add_argument(
        "-mod", "--merge-overlapping-domain",
        type=int,
        choices=[0, 1],
        default=1,
        help="Specify whether to merge overlapping domain hits (1 = merge, 0 = do not merge). Default is 1."
    )

    parser.add_argument(
        "-i", "--input_file",
        type=str,
        required=True,
        help="Required: Path to the input CDD file."
    )

    parser.add_argument(
        "-d", "--debug",
        type=int,
        choices=[0, 1],
        default=0,
        help="Enable debug logs (1 = enable, 0 = disable). Default is 0."
    )

    parser.add_argument(
        "-t", "--tcids",
        type=str,
        default="",
        help="Take in a targeted list of TCIDs to process. \n\nInput can be a file (one family TCID per line), or comma-seperated list (i.e. 1.A.12,1.C.105). Default is all families in cdd file."
    )

    args = parser.parse_args()
    isFile = True
    # Validate input file
    if not os.path.isfile(args.input_file):
        parser.error(f"The input file '{args.input_file}' does not exist.")
    
    # Check Target isFile
    if not os.path.isfile(args.tcids):
        isFile = False

    return args, isFile

def main():
    """
    Main execution function that processes protein domain data.
    
    Workflow:
    1. Parses command-line arguments
    2. Reads and cleans CDD file data
    3. Processes target TCIDs (either from file, comma-separated list, or all)
    4. Validates family IDs
    5. Processes each family to generate domain architecture analysis
    6. Outputs results to output.csv
    """
    # Extract argument values
    args, isfile = parse_arguments()

    # Configuration settings from arguments
    merge = bool(args.merge_overlapping_domain)
    input_file = args.input_file
    debug = bool(args.debug)
    target_ids = args.tcids
    families = []

    # Parse CDD File and extract unique family IDs
    clean = get_clean(input_file)
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    valid_famids = clean["family"].unique()

    # Process target TCIDs based on input method
    if args.tcids:
        if isfile:
            # Read families from file (one per line)
            with open(target_ids, "r") as file:
                for line in file:
                    families.append(line.strip())
        else:
            # Parse comma-separated list of TCIDs
            families = [tcid.strip() for tcid in target_ids.split(",")]
    else:
        # Use all families from CDD file if no targets specified
        families = valid_famids

    # Remove duplicates and sort families
    families = sorted(families)
    invalids = set()

    # Validate family IDs format and existence
    for famid in families:
        famid_split = famid.split(".")
        if len(famid_split) != 3:
            print(f"{famid} is not a valid family TCID.")
            continue

        if famid not in valid_famids:
            print(f"{famid} is not found, skipping...")
            invalids.add(famid)

    # Track execution time
    import time
    start = time.time()
    
    # Remove invalid families and sort
    families = set(families) - invalids
    families = sorted(list(families))

    # Use all valid families if none were selected
    if len(families) == 0:
        families = valid_famids

    # Process each family and generate output rows
    rows = []
    print('\n\n\nProcessing Families...\n\n')
    for cnt, fam in enumerate(families):
        curr_fam = Family(clean[clean["family"] == fam], fam)
        # curr_fam.plot_char()
        # curr_fam.plot_general()
        # curr_fam.plot_summary()
        # curr_fam.plot_holes()
        # curr_fam.plot_arch()
        # curr_fam.generate_checklist()
        for row in curr_fam.generate_csv_rows():
            rows.append(row)
        # Display progress
        print("Processing", fam, str(round(float((cnt + 1) * 100) / len(families), 2)) + "%")

    # Write results to CSV file
    with open("output.csv", "w") as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(["Accession", "Length", "Family", "Subfamily", "Domains", "Seperators"])
        writer.writerows(rows)

    # Display execution time
    end = time.time()
    print("Process took", int((end - start) // 60), "minutes", int(round(end - start, 0)) % 60, "seconds.")

if __name__ == '__main__':
    main()