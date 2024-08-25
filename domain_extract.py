'''
 * Filename: domain_extract.py
 * Author: Xiangyu(Leo) Shi
 * Created: Aug 2024
'''

from util import *
import argparse
import numpy as np

def system_query(clean, OUT_FILE, STD, domains, max_domain):
    accessions = clean['query acc.'].unique()
    #Filter characteristic domains
    dom_df = clean[clean['subject accs.'].isin(domains)]
    cnt = 0
    for acc in accessions:
        print(cnt, acc)

        #Find target accession
        acc_df = dom_df[dom_df['query acc.'] == acc]

        with open(OUT_FILE, 'a') as file:
            # Case 1: No hit found
            if acc_df.shape[0] == 0:
            
                file.write("!" + padding(acc.split('-')[0]))
                file.write('\u2015' * SCALE + "\n")
                continue

            
            acc_len = acc_df['query length'].unique()[0]
            # Get gaps
            domain_regions = np.array(acc_df[['q. start', 'q. end']].sort_values('q. start'))
            links = find_links(acc_len, domain_regions)

            plot = None
            max_int = acc_df[acc_df['subject accs.'] == max_domain][['q. start', 'q. end']].sort_values('q. start')

            # Case 2: No max characteristic domain found
            if max_int.shape[0] == 0:
                plot = plot_structure(acc_len, links)
                plot = plot_scale(plot, STD)
                print(acc + " no char domain.")
                file.write("!" + padding(acc.split('-')[0]))
            else:
                max_int = np.array(max_int)[0]
                plot = plot_structure(acc_len, links)
                plot = plot_special(plot, max_int)
                plot = plot_scale(plot, STD)
                file.write(">" + padding(acc.split('-')[0]))
            #file.write("Unmapped links:\n" + str(links) + "\n")
            file.write(plot + "\n")
        cnt += 1

def family_query(clean, OUT_FILE, STD):
    clean['family'] = clean['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))

    #Filter
    #clean = clean[clean['family'] == '1.B.40']

    families = clean['family'].unique()
    for fam in families:
        fam_df = clean[clean['family'] == fam]
        char_dms, max_dm = char_domains(fam_df)
        print(max_dm)
        with open(OUT_FILE, 'a') as file:
            file.write("\n\n")
            file.write("#" + fam + "\n")
            if len(char_dms) > 0:
                file.write("Characteristic Domains:\n" + str(char_dms) + "\n")
            else:
                file.write("POTENTIAL CLASSIFICATION ERROR! \n")
        file.close()
        system_query(fam_df, OUT_FILE, STD, char_dms, max_dm)



def main():
    parser = argparse.ArgumentParser(description="Extracts domains from given .cdd file.")
    parser.add_argument("input_file", type=str, help="Input file path")
    parser.add_argument("output_file", type=str, help="Output file path")
    parser.add_argument("-s", "--standardize", type=int, default=1, help="Structure plot standardizing, true(1) by default.")

    args = parser.parse_args()
    IN_FILE = args.input_file
    OUT_FILE = args.output_file
    STD = bool(args.standardize)

    clean = get_clean(IN_FILE)
    clean['subject accs.'] = clean['subject accs.'].apply(lambda x: x[4:])

    with open(OUT_FILE, 'w') as file:
        file.write("")
    file.close()

    family_query(clean, OUT_FILE, STD)

    file.close()
if __name__ == '__main__':
    main()