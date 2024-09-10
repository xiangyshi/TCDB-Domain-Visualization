'''
 * Filename: domain_extract.py
 * Author: Xiangyu(Leo) Shi
 * Created: Aug 2024
'''

from util import *
from config import *
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Extracts domains from given .cdd file.")
    parser.add_argument("-m", "--merge", type=int, default=0, help="Merge overlapping domain hits.")

    args = parser.parse_args()
    MERGE = bool(args.merge)

    import time

    # Calculate the start time
    start = time.time()


    clean = get_clean(INPUT_FILE) #TODO Preferablly merge here
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    #clean = clean[clean['query acc.'] == '1.B.40.1.1-P0C2W0']
    with open(ERROR_FILE, 'w') as file:
        file.write('')
        
    families = clean['family'].unique()
    tot_fams = len(families)
    with open('error.txt', 'r') as file:
        for line in file:
            if line[0] != '>':
                continue
            fam = line[1:].strip()
            # print('Generating plots for', fam, str(round(i * 100 / tot_fams, 2))+ "%")
            print(fam)
            curr_fam_data = clean[clean['family'] == fam]
            test_fam = Family(curr_fam_data, fam , MERGE)
            test_fam.plot_general()

    end = time.time()
    print("Plot generation took:", round(end - start, 2), "seconds.")

    #test_fam.plot_general2(palette)


    
    # clean['subject accs.'] = clean['subject accs.'].apply(lambda x: x[4:])

    # with open(OUT_FILE, 'w') as file:
    #     file.write("")
    # file.close()

    # family_query(clean, OUT_FILE, STD)

    # file.close()
if __name__ == '__main__':
    main()