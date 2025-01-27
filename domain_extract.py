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

    # Extract argument values
    args, isfile = parse_arguments()

    merge = bool(args.merge_overlapping_domain)
    input_file = args.input_file
    debug = bool(args.debug)
    target_ids = args.tcids
    families = []

    # Parse CDD File
    clean = get_clean(input_file)
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    valid_famids = clean["family"].unique()

    if args.tcids:
        # Parse Targeted TCIDs
        if isfile:
            # Parse families file
            with open(target_ids, "r") as file:
                for line in file:
                    families.append(line.strip())
        else:
            families = [tcid.strip() for tcid in target_ids.split(",")]
    else:
        families = valid_famids

    # Remove Duplicates
    families = sorted(families)
    invalids = set()

    # Validate IDs
    for famid in families:
        famid_split = famid.split(".")
        if len(famid_split) != 3:
            print(f"{famid} is not a valid family TCID.")
            continue

        if famid not in valid_famids:
            print(f"{famid} is not found, skipping...")
            invalids.add(famid)

    # Main
    import time

    # Calculate the start time
    start = time.time()
    
    families = set(families) - invalids
    families = sorted(list(families))

    # Check No Select
    if len(families) == 0:
        families = valid_famids

    # Process Families
    print('\n\n\nProcessing Families...\n\n')
    for cnt, fam in enumerate(families):
        curr_fam = Family(clean[clean["family"] == fam], fam)
        curr_fam.plot_char()
        curr_fam.plot_general()
        curr_fam.plot_summary()
        curr_fam.plot_holes()
        curr_fam.plot_arch()
        curr_fam.generate_checklist()
        print("Processing", fam, str(round(float((cnt + 1) * 100) / len(families), 2)) + "%")

    end = time.time()
    print("Process took", int((end - start) // 60), "minutes", int(round(end - start, 0)) % 60, "seconds.")


if __name__ == '__main__':
    main()