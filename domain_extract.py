'''
Filename: domain_extract.py
Author: Xiangyu(Leo) Shi
Created: Aug 2024

This script processes protein domains using rpsblast and TCDB-specific proteins.
It analyzes domain data from CDD files and generates domain architecture analysis
for specified protein families.
'''

from utility.util import *
from utility.cdd_parser import clean_cdd
from utility.rescue_parser import clean_rescue
from utility.config import *
import argparse
import numpy as np
import pandas as pd
import csv
from Family import *
from RescueFamily import *

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
        "-m", "--merge-overlapping-domain",
        type=int,
        choices=[0, 1],
        default=1,
        help="Specify whether to merge overlapping domain hits (1 = merge, 0 = do not merge). Default is 1."
    )

    # Create a mutually exclusive group for input types
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-c", "--cdd_input",
        type=str,
        help="Path to the input CDD file. Cannot be used with --rescue_input."
    )
    input_group.add_argument(
        "-r", "--rescue_input",
        type=str,
        help="Path to folder containing rescue files. Cannot be used with --cdd_input."
    )

    parser.add_argument(
        "-p", "--plot",
        type=str,
        default=None,
        help="Plots to generate. CDD options: all / general, char, arch, holes, summary. Rescue options: all / char_rescue. Default: None"
    )

    parser.add_argument(
        "-d", "--data",
        type=str,
        default="",
        help="Output csv file to save protein architecture data for AI Classification. Default: None"
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        default="output/",
        help="Output directory. Default: output/"
    )

    parser.add_argument(
        "-t", "--tcids",
        type=str,
        default="",
        help="Take in a targeted list of TCIDs to process. \n\nInput can be a file (one family TCID per line), or comma-seperated list (i.e. 1.A.12,1.C.105). If not provided, all families in cdd file / rescue folder will be processed."
    )

    args = parser.parse_args()
    isFile = True
    
    # Validate input file/directory
    if args.cdd_input and not os.path.isfile(args.cdd_input):
        parser.error(f"The CDD input file '{args.cdd_input}' does not exist.")
    if args.rescue_input and not os.path.isdir(args.rescue_input):
        parser.error(f"The rescue input path '{args.rescue_input}' is not a directory.")
    
    # Check Target isFile
    if not os.path.isfile(args.tcids):
        print(f"The argument -t/--tcids is not a file: {args.tcids}, treating as comma-separated list.")
        isFile = False

    return args, isFile

def main():
    """
    Main execution function that processes protein domain data.
    
    Workflow:
    1. Parses command-line arguments
    2. Reads and cleans CDD file data or processes rescue IDs
    3. Processes target TCIDs (either from file, comma-separated list, or all)
    4. Validates family IDs
    5. Processes each family to generate domain architecture analysis
    6. Outputs results to output.csv
    """
    # Extract argument values
    args, isfile = parse_arguments()

    # Configuration settings from arguments
    merge = bool(args.merge_overlapping_domain)
    plot_options = args.plot.split(",") if args.plot is not None else []
    datafile = args.data if args.data else None
    outputdir = args.output

    # Create output directory if it doesn't exist
    os.makedirs(outputdir, exist_ok=True)

    process_fam_ids = set()

    # Process target TCIDs based on input method
    if args.tcids:
        if isfile:
            # Read families from file (one per line)
            with open(args.tcid, "r") as file:
                for line in file:
                    process_fam_ids.add(line.strip())
        else:
            # Parse comma-separated list of TCIDs
            process_fam_ids = [tcid.strip() for tcid in args.tcids.split(",")]

    # Parse input based on which input was provided
    if args.cdd_input:
        # CDD file mode
        family_data = clean_cdd(args.cdd_input, process_fam_ids)
    else:
        family_data = clean_rescue(args.rescue_input, process_fam_ids)

    # Track execution time
    import time
    start = time.time()
    
    unique_fam_ids = family_data['family'].unique()

    # Process each family and generate output rows
    rows = []
    print('\n\n\nProcessing Families...\n\n')

    for cnt, fam in enumerate(unique_fam_ids):
        print("Processing", fam, str(round(float((cnt) * 100) / len(unique_fam_ids), 2)) + "%")
        if args.rescue_input:
            curr_fam = RescueFamily(family_data[family_data['family'] == fam], fam, outputdir, merge)
            curr_fam.plot_char_rescue()
        else:
            curr_fam = Family(family_data[family_data['family'] == fam], fam, outputdir, merge)

        if "all" in plot_options:
            curr_fam.plot()
        else:
            curr_fam.plot(options=plot_options)

        if args.data:
            curr_data = curr_fam.generate_csv_rows()
            rows += curr_data

    
    # Write results to CSV file
    with open(os.path.join(outputdir, datafile), "w") as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(["Accession", "Length", "Family", "Subfamily", "Domains", "Seperators"])
        writer.writerows(rows)

    # Display execution time
    end = time.time()
    print("Process took", int((end - start) // 60), "minutes", int(round(end - start, 0)) % 60, "seconds.")

if __name__ == '__main__':
    main()