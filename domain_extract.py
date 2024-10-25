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
    parser.add_argument("-m", "--merge", type=int, default=1, help="Merge overlapping domain hits.")
    parser.add_argument("-i", "--input_file", type=str, help="Path of input CDD file.")
    parser.add_argument("-d", "--debug", type=int, default=0, help="Show debug logs.")

    args = parser.parse_args()
    MERGE = bool(args.merge)
    INPUT_FILE = args.input_file
    DEBUG = bool(args.debug)

    # Main
    import time

    # Calculate the start time
    start = time.time()


    clean = get_clean(INPUT_FILE) #TODO Preferablly merge here
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    # clean = clean[clean['family'] == "1.A.119"]
    # clean = clean[clean['query acc.'] == '1.A.104.1.1-P76298']
    # with open(ERROR_FILE, 'w') as file:
    #     file.write('')

         
    families = clean['family'].unique()
    for cnt, fam in enumerate(families):
        curr_fam = Family(clean[clean["family"] == fam], fam)
        curr_fam.plot_arch()
        curr_fam.plot_char()
        curr_fam.plot_general()
        curr_fam.plot_summary()
        curr_fam.plot_holes()
        print(fam, str(round(float(cnt * 100) / len(families), 2)) + "%")
        #test_fam = Family(clean[clean['family'] == fam], fam)
        #test_fam.plot_holes()
    # end = time.time()
    # print("Process took", int((end - start) // 60), "minutes", int(round(end - start, 0)) % 60, "seconds.")

    #test_fam.plot_general2(palette)


    
    # clean['subject accs.'] = clean['subject accs.'].apply(lambda x: x[4:])

    # with open(OUT_FILE, 'w') as file:
    #     file.write("")
    # file.close()

    # family_query(clean, OUT_FILE, STD)

    # file.close()
if __name__ == '__main__':
    main()