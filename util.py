# Utility
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

# Plotting
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import mpld3

class Family:
    """
    The Family class is used primarily to facilitate the organization of TCDB
    and plotting systems collectively with standardization.

    Attributes:
        data (DataFrame): Contains the cleaned output of this family from rpblast.
        fam_id (String): Unique TCDB accession of this family.
        systems (System[]): A collection of systems in the family represented by the System class.
        char_domains (String[]): The set of characteristic domains of this family.
        max_sys_len (Integer): The length of the longest protein in the family (for plot standardization).

    Methods:
        Family(data: DataFrame, fam_id: String): Instantiates a Family object.
        get_systems(): Returns the list of systems of this family.
        get_domains(): TODO
        get_char_domains(): Returns the list of characteristic domains of this family.
        plot_general(): Generates the general plot in html format.
    """
     
    def __init__(self, data, fam_id, merge=True):
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

        with open(ERROR_FILE, 'a') as file:
            file.write(self.fam_id + '  ' + str(bad_sys) + '\n')

    
    def get_systems(self):
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
        domains = data[[DOM_START, DOM_END, DOM_ID]].sort_values(DOM_START)
        domains_loc = np.array(domains[[DOM_START, DOM_END]]).tolist()
        links = find_links(sys_len, domains_loc)

        domains = links + np.array(domains).tolist()
        domains = sorted(domains, key=lambda x: x[0])
        return domains
    
    def get_char_domains(self, thresh=0.5):
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

        gen_doms = set(self.data[DOM_ID].unique()) - set(self.char_domains)
        gen_palette = sns.color_palette("husl", n_colors=len(gen_doms))
        gen_palette = {domain: mcolors.to_hex(gen_palette[i]) for i, domain in enumerate(gen_doms)}

        char_palette = sns.color_palette("husl", n_colors=len(self.char_domains))
        char_palette = {domain: mcolors.to_hex(char_palette[i]) for i, domain in enumerate(self.char_domains)}

        return gen_palette, char_palette
    
    def plot_palette(self):
        # Get the list of colors
        colors = list(self.gen_palette.values())
        
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(10, 2))
        
        # Create a color bar
        n = len(colors)
        for i, color in enumerate(colors):
            ax.add_patch(plt.Rectangle((i / n, 0), 1 / n, 1, color=color))
        
        # Set axis limits and labels
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Color Palette')
        
        # Display the plot
        plt.show()

    def plot_general(self, mode='char'):
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
        combine_svgs(svgs, 'plot_character.html')

    def plot_summary(self):
        fig, ax = plt.subplots(figsize=(20, 9))
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

class System(Family):
    def __init__(self, fam_id, sys_id, sys_len, domains):
        self.fam_id = fam_id
        self.sys_id = sys_id
        self.sys_len = sys_len
        self.domains = domains
    
    def check_char(self, char_domains):
        for dom in self.domains:
            if dom[2] in char_domains:
                self.has_char = True
                return
        self.has_char = False

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