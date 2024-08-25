import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from config import CLEAN_COMMAND, CUT_COMMAND, SCALE

def padding(s):
        res = s
        while len(res) < 16:
            res += ' '
        return res

def get_clean(in_file) -> pd.DataFrame:
    '''
    Cleans the input data by removing comments and selecting specific fields,
    then returns the cleaned data as a pandas DataFrame.
    '''
    try:
        # Obtain column headers
        header = ""
        with open(in_file, mode='r') as input_file:
            for line in input_file:
                if line.startswith("# Fields:"):
                    header = line.strip().replace("# Fields: ", "")  # Extract and clean the header line
                    break  # Exit after finding the header

        if not header:
            raise ValueError("Fields header not found in the input file.")
        
        # Run the CLEAN command
        clean_process = subprocess.run(CLEAN_COMMAND, capture_output=True, text=True)

        if clean_process.returncode != 0:
            raise Exception(f"Error in CLEAN command: {clean_process.stderr}")

        # Run the CUT command on the output of the CLEAN command
        cut_process = subprocess.run(CUT_COMMAND, input=clean_process.stdout, capture_output=True, text=True)

        if cut_process.returncode != 0:
            raise Exception(f"Error in CUT command: {cut_process.stderr}")

        # Convert output to DataFrame
        from io import StringIO
        cleaned_data = StringIO(cut_process.stdout)
        df = pd.read_csv(cleaned_data, sep='\t', header=None)

        # Set column names
        header = header.split(", ")
        header = header[0:5] + header[6:8]
        df.columns = header

        return df  # Return the final output
    
    except Exception as e:
        print("An error occurred:", e)
        return None

def find_links(p_length, domain_regions):
    unknown = [1] * p_length
    
    for start, end in domain_regions:
        for i in range(start - 1, end):
            #try:
            unknown[i] = 0
            #except Exception as e:
            #    print(i, p_length, start, end)
    res = []
    start = None

    for i in range(p_length):
        if unknown[i] == 1:
            if start is None:
                start = i + 1  # Convert to 1-based index
        else:
            if start is not None:
                res.append([start, i])  # Use i+1 to account for 1-based index
                start = None

    if start is not None:
        res.append([start, p_length])
    return res

def char_domains(fam_df, thresh=0.5):
    systems = fam_df['query acc.'].unique()
    domains = fam_df['subject accs.'].unique()
    domain_counts = []
    for system in systems:
        sys_df = fam_df[fam_df['query acc.'] == system]
        for domain in domains:
            if np.any(domain == sys_df['subject accs.']):
                domain_counts.append(domain)
    domain_counts = Counter(domain_counts)
    res = []
    max_domain = ""
    max_cnt = 0
    for dom, cnt in domain_counts.items():
        if cnt / len(systems) > thresh:
            res.append(dom)
            if max_cnt < cnt:
                max_cnt = cnt
                max_domain = dom
    return res, max_domain

def plot_structure(p_length, links):
    structure = []
    curr = 1
    while len(links) != 0:
        if links[0][0] == curr:
            structure.append([*links.pop(0), 0])
            curr = structure[-1][1] + 1
        else:
            structure.append([curr, links[0][0] - 1, 1])
            curr = links[0][0]
    if len(structure) == 0:
        structure.append([0, p_length, 1])
    elif structure[-1][1] != p_length:
        if structure[-1][2] == 1:
            raise Exception("Structure Plot Exception")
        structure.append([structure[-1][1] + 1, p_length, 1])

    res = ''
    for curr in structure:
        while len(res) < curr[1]:
            if curr[-1] == 0:
                res += '\u2015'
            else:
                res += '≡'
    return res

def plot_special(structure, max_domain=[0, 0]):
    res = ""
    if np.any(max_domain):
        start = max_domain[0]
        end = max_domain[1]
        for i in range(len(structure)):
            if i < start or i > end:
                res += structure[i]
            else:
                res += '\u2588'
    return res

def plot_scale(structure, scaled=True):
    scale = SCALE
    if not scaled:
        return structure
    res = ''
    for n in range(scale):
        res += structure[int(np.floor(n * len(structure) / scale))]
    return res

def plot_domains(intervals):
    for i in range(len(intervals)):
        space = np.linspace(intervals[i][0], intervals[i][1], 2)
        plt.plot(space, [i] * len(space))
    plt.show()
