'''
Module Name: util.py

Description:
    This module provides utility functions and classes to the driver 
    code `domain_extract.py`, including the following classes:

        `Family`: represents a protein family.
        `System`: represents a system of a protein family.
            (currently only supports single component systems)
        `Domain`: represents a domain in a sequence of protein.
'''
import subprocess
import platform
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.transforms as mtrs
import seaborn as sns
from collections import Counter
from config import *
import scipy.stats as stats
import networkx as nx

# Plotting
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

class Family:
    '''
    The Family class is used primarily to facilitate the organization of TCDB
    and plotting systems collectively with standardization.

   
    Methods:
        Family(data: DataFrame, fam_id: String): Instantiates a Family object.
        get_systems(): Returns the list of systems of this family.
        get_domains(): TODO
        get_char_domains(): Returns the list of characteristic domains of this family.
        plot_general(): Generates the general plot in html format.
    '''
     
    def __init__(self, data, fam_id, merge=True):
        '''
            Extracts and constructs a Family object.

            Attributes:
            data (pd.DataFrame): Contains the cleaned output of this family from rpblast.
            fam_id (String): Unique TCDB accession of this family.
            systems (System (list)): A collection of systems in the family represented by the System class.
            char_domains (String (list)): The set of characteristic domains of this family.
            max_sys_len (Integer): The length of the longest protein in the family (for plot standardization).
        '''
        self.data = data
        self.fam_id = fam_id
        self.merge = merge
        #print("Gathering systems.")
        self.systems = self.get_systems()
        #print("Gathering char domains.")
        self.char_domains = self.get_char_domains()
        #print("Checking systems without char_domains.")
        bad_sys = []
        for sys in self.systems:
            sys.check_char(self.char_domains)
            if not sys.has_char:
                bad_sys.append(sys.sys_id)
        #print("Creating palette.")
        self.gen_palette, self.char_palette = self.get_palettes()
    
    def get_systems(self):
        '''
        Retrieves the different systems stored in `self.data`, 
        
        Returns:
            System (list): a list of constructed System objects.
        '''
        res = []
        sys_ids = self.data[SYS_ID].unique()
        self.max_sys_len = 0
        for sys_id in sys_ids: 
            curr = self.data[self.data[SYS_ID] == sys_id]
            sys_len = curr[SYS_LEN].unique()[0]
            if sys_len > self.max_sys_len:
                self.max_sys_len = sys_len
            domains = self.get_domains(sys_len, curr)
            if self.merge:
                domains = merge_domains(domains)
            res.append(System(self.fam_id, sys_id, sys_len, domains))
        return res
    
    def get_domains(self, sys_len, data):
        '''
        Retrieves the list of domains in a given system. Requires 
        data to contain domain hits from the same system.
        
        Attributes:
            sys_len (int): the total length of the system.
            data (pd.DataFrame): data containing hits only within this system.

        Returns:
            list: A list representation of domain hits, including
                links that are present.
        '''
        domains = data[[DOM_START, DOM_END, DOM_ID]].sort_values(DOM_START)
        domains_loc = np.array(domains[[DOM_START, DOM_END]]).tolist()
        links = find_links(sys_len, domains_loc)

        domains = links + np.array(domains).tolist()
        domains = sorted(domains, key=lambda x: x[0])
        return domains
    
    def get_char_domains(self, thresh=0.5):
        '''
        Retrieves a list of characteristic domains of the given family.

        Attributes:
            thresh (float): A value between 0 and 1, represents the threshold
                precence percentage for a domain to be considered characteristic.
                (i.e. 0.5 means if a domain is present in more than 50% of the
                systems in the family, then it is characteristic.)
        
        Returns:
            list: A list of domain id's that are considered characteristic.
        '''
        systems = self.data[SYS_ID].unique()
        domains = self.data[DOM_ID].unique()
        domain_counts = []
        for system in systems:
            sys_df = self.data[self.data[SYS_ID] == system]
            for domain in domains:
                if np.any(domain == sys_df[DOM_ID]):
                    domain_counts.append(domain)
        domain_counts = Counter(domain_counts)
        res = []
        for dom, cnt in domain_counts.items():
            if cnt / len(systems) > thresh:
                res.append(dom)
        return res

    def get_palettes(self):
        '''
        A helper function to create a color palette for enhanced visualization.
        '''
        gen_doms = set(self.data[DOM_ID].unique()) - set(self.char_domains)
        gen_palette = sns.color_palette("husl", n_colors=len(gen_doms))
        gen_palette = {domain: mcolors.to_hex(gen_palette[i]) for i, domain in enumerate(gen_doms)}

        char_palette = sns.color_palette("husl", n_colors=len(self.char_domains))
        char_palette = {domain: mcolors.to_hex(char_palette[i]) for i, domain in enumerate(self.char_domains)}

        return gen_palette, char_palette

    def plot_general(self, mode='char'):
        '''
        '''
        svgs = []
        for i, sys in enumerate(self.systems):
            size = len([0 for i in sys.domains if i[-1] != "-1"])
            fig, ax = plt.subplots(figsize=(16, 0.25 * (size + 2)))  # Adjust size as needed
            
            ax.set_title(sys.sys_id)
            ax.set_xlabel('Residual')
            ax.set_ylabel('Domains')

            # Plot domains
            space = np.linspace(0, sys.sys_len - 1, 2)
            cnt = 0
            dom_ids = [sys.sys_id.split('-')[0]]
            ax.plot(space, [cnt] * 2)

            # Set plotting mode (char first or mixed)
            gen_doms = sys.domains
            if mode == 'char':
                gen_doms = [sys for sys in sys.domains if sys[-1] in self.char_domains] + [sys for sys in sys.domains if sys[-1] not in self.char_domains]

            for dom in gen_doms:
                if dom[-1] == "-1":
                    continue
                cnt -= 1
                dom_ids.append(dom[2])
                space = np.linspace(dom[0], dom[1], 2)
                if dom[2] in self.gen_palette:
                    ax.plot(space, [cnt] * 2, color=self.gen_palette[dom[2]], label=dom[2], linewidth=8)
                else:
                    ax.plot(space, [cnt] * 2, color=CHAR_COLOR, label=dom[2], linewidth=8)

            # Set consistent axis limits
            ax.set_xlim(0, self.max_sys_len)  # Adjust based on your data
            ax.set_ylim(cnt - 1, 1)
            
            ax.set_yticks(range(0, cnt - 1, -1))
            ax.set_yticklabels(dom_ids)

            # Save as SVG
            svg_file = f"plots/plot_{sys.sys_id}.svg"
            fig.savefig(svg_file, format='svg', bbox_inches='tight', pad_inches=0.2)
            svgs.append(svg_file)
            plt.close(fig)  # Close the figure to avoid memory issues
        combine_svgs(svgs, "general/"  + self.fam_id + "-e4.html")
    
    def plot_char(self):
        svgs = []
        for i, sys in enumerate(self.systems):
            char_doms = [dom for dom in sys.domains if dom[-1] in self.char_domains]
            fig, ax = plt.subplots(figsize=(16, 0.25 * (len(char_doms) + 2)))  # Adjust size as needed
            
            ax.set_title(sys.sys_id)
            ax.set_xlabel('Residual')
            ax.set_ylabel('Domains')

            # Plot domains
            space = np.linspace(0, sys.sys_len - 1, 2)
            cnt = 0
            dom_ids = [sys.sys_id.split('-')[0]]
            ax.plot(space, [cnt] * 2)

            for dom in char_doms:
                cnt -= 1
                dom_ids.append(dom[2])
                space = np.linspace(dom[0], dom[1], 2)
                ax.plot(space, [cnt] * 2, color=self.char_palette[dom[2]], label=dom[2], linewidth=8)

            # Set consistent axis limits
            ax.set_xlim(0, self.max_sys_len)  # Adjust based on your data
            ax.set_ylim(cnt - 1, 1)
            
            ax.set_yticks(range(0, cnt - 1, -1))
            ax.set_yticklabels(dom_ids)

            # Save as SVG
            svg_file = f"plots/plot_{sys.sys_id}.svg"
            fig.savefig(svg_file, format='svg', bbox_inches='tight', pad_inches=0.2)
            svgs.append(svg_file)

            plt.close(fig)  # Close the figure to avoid memory issues
        combine_svgs(svgs, '/char/' + self.fam_id + "-char.html")

    def plot_summary(self):
        fig, ax = plt.subplots(figsize=(16, 0.25 * (len(self.systems) + 2)))
        for i, sys in enumerate(self.systems):
            for dom in sys.domains:
                if dom[-1] == "-1":
                    ax.barh(i, dom[1] - dom[0], left=dom[0], height=0.03, color='blue')
                else:
                    ax.barh(i, dom[1] - dom[0], left=dom[0], height=0.4, color='blue')

        # Set y-ticks and y-tick labels
        ax.set_yticks(range(len(self.systems)))  # Ensure y-ticks match the number of systems
        sys_ids = [sys.sys_id.split('-')[0] for sys in self.systems]  # Collect the system IDs
        ax.set_yticklabels(sys_ids)  # Set the custom y-tick labels
        
        # Set axis labels and title
        ax.set_xlabel('Residual')
        ax.set_title(self.fam_id + ' Summary')

        # Save the plot to an HTML file
        fig.savefig("plots/summary/"  + self.fam_id + "-summary.svg")
        plt.close()

    def plot_arch(self):
        if len(self.systems) <= 2:
            print("System count too small, unable to accurately construct architecture. Skipping arch plot...")
            return
        
        char_domain_locs = {k:[] for k in self.char_domains}
        char_domain_map = {}
        # Store char domain locations in format "char_dom_id" : [[#char_dom_id appearances in sys1], [#char_dom_id appearances in sys2], ... , [#char_dom_id appearances in sysn]]
        for system in self.systems:
            sys_dom_locs = {k:[] for k in self.char_domains}
            for domain in system.domains:
                if domain[-1] in char_domain_locs:
                    sys_dom_locs[domain[-1]].append([domain[0] / system.sys_len, domain[1] / system.sys_len])
            
            
            for char_dom in self.char_domains:
                char_domain_locs[char_dom].append(sys_dom_locs[char_dom])
        # Analyze structure, choose majority appearances (mode), update
        for char_dom in self.char_domains:
            curr_lst = char_domain_locs[char_dom]
            
            # TODO: better detection than mode
            curr_mode = stats.mode([len(lst) for lst in curr_lst if len(lst) > 0])
            curr_lst = [lst for lst in curr_lst if len(lst) == curr_mode.mode]
            if curr_mode.mode > 1:
                for i in range(curr_mode.mode):
                    char_domain_locs[char_dom + ' ' + str(i)] = [lst[i] for lst in curr_lst]
                del char_domain_locs[char_dom]
            else:
                char_domain_locs[char_dom] = [lst[0] for lst in curr_lst]

        for char_dom, intervals in char_domain_locs.items():
            start = confidence_interval_mean([i[0] for i in intervals])
            end = confidence_interval_mean([i[1] for i in intervals])
            
            components = char_dom.split(' ')
            if len(components) > 1:
                if components[0] in char_domain_map:
                    char_domain_map[components[0]].append([start * 100, end * 100])
                else:
                    char_domain_map[components[0]] = [[start * 100, end * 100]]
            else:
                char_domain_map[char_dom] = [[start * 100, end * 100]]
    
        #plot
        cnt = 0
        fig, ax = plt.subplots(figsize=(16, 0.25 * (len(self.char_domains) + 2)))
        for char_dom, intervals in char_domain_map.items():
            for interval in intervals:
                ax.barh(cnt, interval[1]-interval[0], left=interval[0], height=0.4, color=self.char_palette[char_dom.split(' ')[0]])
            cnt += 1
        
        ax.set_yticks(range(len(char_domain_map)))  # Ensure y-ticks match the number of systems
        ax.set_yticklabels(list(char_domain_map.keys()))  # Set the custom y-tick labels
        
        ax.set_xlabel('Residual Percentile %')
        ax.set_title(self.fam_id + ' Architecture')

        ax.set_xlim(-1, 101)

        svg = "plots/arch/"  + self.fam_id + "-arch.svg"
        fig.savefig(svg)
        plt.close()
        combine_svgs([svg], '/arch/' + self.fam_id + "-arch.html")
    
    def plot_holes(self, thresh=50):
        fig, ax = plt.subplots(figsize=(16, 0.25 * (len(self.systems) + 2)))
        for i, sys in enumerate(self.systems):
            ax.barh(i, sys.sys_len, left=0, height=0.03, color='red')
            for dom in sys.domains:
                if dom[-1] == "-1" and dom[1] - dom[0] >= thresh:
                    ax.barh(i, dom[1] - dom[0], left=dom[0], height=0.4, color='red')
            

        # Set y-ticks and y-tick labels
        ax.set_yticks(range(len(self.systems)))  # Ensure y-ticks match the number of systems
        sys_ids = [sys.sys_id.split('-')[0] for sys in self.systems]  # Collect the system IDs
        ax.set_yticklabels(sys_ids)  # Set the custom y-tick labels
        
        # Set axis labels and title
        ax.set_xlabel('Residual')
        ax.set_title(self.fam_id + ' Holes')

        # Save the plot to an HTML file
        fig.savefig("plots/holes/"  + self.fam_id + "-holes.svg")
        plt.close()
    
    def get_holes(self):
        holes_dict = {sys.sys_id:sys.holes for sys in self.systems}
        if sum([len(val) for val in holes_dict.values()]) == 0:
            print(f"Family {self.fam_id} have no holes.")
            return
        return holes_dict

    def generate_checklist(self):
        holes = []
        pairs = []
        for sys in self.systems:
            for hole in sys.holes:
                holes.append(hole)
        print(holes)
        for i in range(len(holes)):
            for j in range(i + 1, len(holes)):
                if compare_margin(holes[i], holes[j]) >= 1:
                    pairs.append([i, j])
        idx_set = get_connected_components(len(holes), pairs)
        
        
        print()
        print()
        print(idx_set)
        res = []
        start = []
        for subset in idx_set:
            res.append([holes[i].sys_id.split('-')[0] for i in subset])
            start.append([holes[i].start for i in subset])
        print(res)
        print(start)
        print()
        print()
        return res
    
class System(Family):
    def __init__(self, fam_id, sys_id, sys_len, domains):
        self.fam_id = fam_id
        self.sys_id = sys_id
        self.sys_len = sys_len
        self.domains = domains
        self.holes = []
        self.get_holes()
    
    def check_char(self, char_domains):
        for dom in self.domains:
            if dom[2] in char_domains:
                self.has_char = True
                return
        self.has_char = False
    
    def get_holes(self, thresh=50, margin=10):
        self.holes = []
        for i, dom in enumerate(self.domains):
            if dom[-1] == "-1" and dom[1] - dom[0] >= thresh:
                left_doms, right_doms = find_margins(self.domains, dom[0] - margin, dom[1] + margin)
                self.holes.append(Hole(self.sys_id, i, dom[0], dom[1], left_doms, right_doms))

class Hole:
    def __init__(self, sys_id, pos, start, end, left, right):
        self.sys_id = sys_id
        self.pos = pos
        self.start = start
        self.end = end
        self.left = left
        self.right = right

class Domain:
    def __init__(self, dom_id, start, end, type):
        self.dom_id = dom_id
        self.start = start
        self.end = end
        self.type = type


def get_clean(in_file) -> pd.DataFrame:
    '''
    Cleans the input data by removing comments and selecting specific fields,
    then returns the cleaned data as a pandas DataFrame.
    '''

    # Adjust command pipeline based on machine architecture
    arch = platform.system()
    CLEAN_COMMAND = UNIX_CLEAN_COMMAND
    if arch == 'Windows':
        CLEAN_COMMAND = WIN_CLEAN_COMMAND
    
    # Append input_file path
    CLEAN_COMMAND += [in_file]

    try:
        # Obtain column headers
        header = None
        with open(in_file, 'r') as input_file:
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

        # Convert output to DataFrame
        from io import StringIO
        cleaned_data = StringIO(clean_process.stdout)
        df = pd.read_csv(cleaned_data, sep='\t', header=None)

        # Set column names
        df.columns = header.split(", ")

        return df
    
    except Exception as e:
        print("An error occurred:", e)
        return None

def find_links(p_length, domain_regions):
    unknown = [1] * p_length
    
    for start, end in domain_regions:
        for i in range(start - 1, end):
            unknown[i] = 0
    res = []
    start = None

    for i in range(p_length):
        if unknown[i] == 1:
            if start is None:
                start = i + 1  # Convert to 1-based index
        else:
            if start is not None:
                res.append([start, i, "-1"])  # Use i+1 to account for 1-based index
                start = None

    if start is not None:
        res.append([start, p_length, "-1"])
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

def combine_svgs(svgs, filename):
        html_content = ''

        for svg in svgs:
            with open(svg, 'r') as file:
                svg_content = file.read().strip()
                # Ensure the SVG content is well-formed and includes the correct XML namespaces
                html_content += f'<div>{svg_content}</div>'
            os.remove(svg)
        html_content += '</body></html>'
        
        with open(f'plots/{filename}', "w") as f:
            f.write(html_content)

def merge_domains(domains):
    # Sort intervals by id and then by the start of the interval
    domains = sorted(domains, key=lambda x: (x[2], x[0]))

    merged = []
    current_start, current_end, current_id = domains[0]

    for i in range(1, len(domains)):
        start, end, id_ = domains[i]

        if id_ == current_id:
            # If the intervals overlap, merge them
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                # If they don't overlap, add the current interval to the result
                merged.append((current_start, current_end, current_id))
                current_start, current_end, current_id = start, end, id_
        else:
            # If the ID changes, add the current interval to the result
            merged.append((current_start, current_end, current_id))
            current_start, current_end, current_id = start, end, id_

    # Add the last interval
    merged.append((current_start, current_end, current_id))
    return sorted(merged, key=lambda x : x[0])

def confidence_interval_mean(coords, confidence_level=0.95):
    if len(set(coords)) == 1:
        return np.mean(coords)
    
    mean = np.mean(coords)
    std_err = stats.sem(coords)
    df = len(coords) - 1

    ci = stats.t.interval(confidence_level, df, loc=mean, scale=std_err)
    return np.mean(ci)

def find_margins(domains, margin_left, margin_right):
    left_doms = []
    right_doms = []
    
    for dom in domains:
        if dom[0] <= margin_left and margin_left <= dom[1]:
            left_doms.append(dom)
        if dom[0] <= margin_right and margin_right <= dom[1]:
            right_doms.append(dom)

    return left_doms, right_doms

def compare_margin(hole1, hole2):
    h1_left = [dom[-1] for dom in hole1.left]
    h2_left = [dom[-1] for dom in hole2.left]
    h1_right = [dom[-1] for dom in hole1.right]
    h2_right = [dom[-1] for dom in hole2.right]
    return len(set(h1_left) & set(h2_left)) + len(set(h1_right) & set(h2_right))

def get_connected_components(n, pair_list):
    G = nx.Graph()
    
    # Add all nodes (1-indexed)
    G.add_nodes_from(range(0, n))
    
    # Add all edges (the pairs from the pair_list)
    G.add_edges_from(pair_list)
    
    # Get connected components
    connected_components = list(nx.connected_components(G))
    
    return connected_components