'''
 * Filename: domain_extract.py
 * Author: Xiangyu(Leo) Shi
 * Created: Aug 2024
'''

from util import *
from config import *
import argparse
import numpy as np
from Family import *
from RescueFamily import *
from rescue_parser import parse_rescue

def parse_arguments():
    """Parses and validates command-line arguments."""
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
        help="Path to the input CDD file."
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
        help="Comma-separated list of TCID families to process (e.g., 1.A.13,1.C.105). Default is 0. If neither -t and -f selected, all proteins will be processed."
    )

    parser.add_argument(
        "-f", "--ftcids",
        type=str,
        default="",
        help="Path to a file containing TCID families to process, one ID per line. If neither -t and -f selected, all proteins will be processed."
    )

    parser.add_argument(
        "-r", "--rescue",
        type=str,
        default="",
        help="Rescue file to parse."
    )

    args = parser.parse_args()

    # Validate input file
    if not os.path.isfile(args.input_file):
        parser.error(f"The input file '{args.input_file}' does not exist.")

    # Validate TCID file if provided
    if args.ftcids and not os.path.isfile(args.ftcids):
        parser.error(f"The TCID file '{args.ftcids}' does not exist.")

    return args

def main():

    # Extract argument values
    args = parse_arguments()

    merge = bool(args.merge_overlapping_domain)
    input_file = args.input_file
    debug = bool(args.debug)
    families = [tcid.strip() for tcid in args.tcids.split(",") if tcid.strip()]
    fam_file = args.ftcids
    resc_file = args.rescue

    # Main
    import time

    # Calculate the start time
    start = time.time()

    # Parse CDD File
    clean = get_clean(input_file)
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    valid_famids = clean["family"].unique()

    # Parse fam_file
    if fam_file:
        with open(fam_file, "r") as file:
            for line in file:
                families.append(line.strip())

    # Remove Duplicates
    families = sorted(families)
    invalids = set()

    # Validate IDs
    for famid in families:
        if famid not in valid_famids:
            print(f"{famid} is not found, skipping...")
            invalids.add(famid)
    
    families = set(families) - invalids
    families = sorted(list(families))

    rescued_domains = parse_rescue(resc_file)
    rescued_domains['family'] = rescued_domains['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    # Check No Select
    if len(families) == 0:
        families = valid_famids

    # Process Families
    print('\n\n\nProcessing Families...\n\n')
    print(clean.columns)
    for cnt, fam in enumerate(families):
        fam_df = clean[clean["family"] == fam]
        resc_df = rescued_domains[rescued_domains["family"] == fam]
        curr_fam = Family(fam_df, fam)
        resc_fam = RescueFamily(resc_df, fam)
        resc_fam.plot_char_rescue()
        # curr_fam.plot_char()
        # curr_fam.plot_general()
        # curr_fam.plot_summary()
        # curr_fam.plot_holes()
        # curr_fam.plot_arch()
        # curr_fam.generate_checklist()
        print("Processing", fam, str(round(float((cnt + 1) * 100) / len(families), 2)) + "%")

    end = time.time()
    print("Process took", int((end - start) // 60), "minutes", int(round(end - start, 0)) % 60, "seconds.")


if __name__ == '__main__':
    main()