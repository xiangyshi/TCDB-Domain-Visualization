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
    parser.add_argument("input_file", type=str, help="Input file path")
    parser.add_argument("output_file", type=str, help="Output file path")
    parser.add_argument("-s", "--standardize", type=int, default=1, help="Structure plot standardizing, true(1) by default.")

    args = parser.parse_args()
    IN_FILE = args.input_file
    OUT_FILE = args.output_file
    STD = bool(args.standardize)
    import time

    # Calculate the start time
    start = time.time()
    clean = get_clean(IN_FILE)
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
    clean = clean[clean['family'] == '1.B.40']
    test_fam = Family(clean, '1.B.40')
    print('Family built.')
    test_fam.plot_general()
    test_fam.plot_char()
    test_fam.plot_summary()
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