import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import seaborn as sns
import mpld3

from collections import Counter
from config import *


class Family:
    def __init__(self, data, fam_id):
        self.data = data
        self.fam_id = fam_id
        print("Gathering systems.")
        self.systems = self.get_systems()
        print("Gathering char domains.")
        self.char_domains = self.get_char_domains()
    
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
            res.append(System(self.fam_id, sys_id, sys_len, domains))
            print('System ', sys_id , sys_len)
        return res
    
    def get_domains(self, sys_len, data):
        domains = data[[DOM_START, DOM_END, DOM_ID]].sort_values(DOM_START)
        domains_loc = np.array(domains[[DOM_START, DOM_END]])
        links = np.array(find_links(sys_len, domains_loc))
        domains = np.array(domains)

        if len(links) == 0:
            links = links.reshape(0, 3)
        if len(domains) == 0:
            domains = domains.reshape(0, 3)

        domains = np.concatenate((links, domains), axis=0)
        return domains[domains[:, 0].argsort()]

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

    def plot_general1(self, palette):
        # Create a figure and axes
        fig, axes = plt.subplots(len(self.systems), 1, figsize=(18, 6 * len(self.systems)))
        if len(self.systems) == 1:
            axes = [axes]

        for i, (sys, ax) in enumerate(zip(self.systems, axes)):
            ax.set_title(sys.sys_id)
            ax.set_xlabel('Residual')
            ax.set_ylabel('Domains')

            # Plot domains
            space = np.linspace(0, sys.sys_len - 1, 2)
            cnt = 0
            ax.plot(space, [cnt] * 2)
            
            for dom in sys.domains:
                if dom[-1] == -1:
                    continue
                cnt -= 1
                space = np.linspace(dom[0], dom[1], 2)
                ax.plot(space, [cnt] * 2, color=palette[dom[2]], label=dom[2], linewidth=8)
            
            # Set consistent axis limits
            ax.set_xlim(0, self.max_sys_len)  # Adjust based on your data
            ax.set_ylim(cnt - 0.5, 0.5)
            #ax.legend()

        # Adjust spacing between plots
        fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.5)
        fig.tight_layout()

        # Save to HTML
        html_str = mpld3.fig_to_html(fig)
        with open("plots/plot1.html", "w") as f:
            f.write(html_str)

    def plot_general2(self, palette):
        # Create a subplot figure
        fig = make_subplots(rows=len(self.systems), cols=1, shared_xaxes=True, vertical_spacing=0.01)

        for i, sys in enumerate(self.systems):
            cnt = 0
            
            # Add the base line
            space = np.linspace(0, sys.sys_len - 1, 2)
            fig.add_trace(go.Scatter(
                x=space,
                y=[cnt] * 2,
                mode='lines',
                name=f'{sys.sys_id} - Base Line',
                line=dict(width=2, color='black')  # Base line color
            ), row=i+1, col=1)

            for dom in sys.domains:
                if dom[-1] == -1:
                    continue
                cnt -= 1
                space = np.linspace(dom[0], dom[1], 2)
                color = palette.get(dom[2], '#000000')  # Get color from palette
                fig.add_trace(go.Scatter(
                    x=space,
                    y=[cnt] * 2,
                    mode='lines',
                    line=dict(color=color, width=6),  # Use color from palette
                    name=dom[2]  # Legend entry
                ), row=i+1, col=1)

                fig.add_annotation(
                    x=0.5,  # Center horizontally
                    y=1.05,  # Position above the plot
                    text=sys.sys_id,  # System ID as title
                    showarrow=False,
                    font=dict(size=14, color="black"),
                    row=i+1,
                    col=1,
                    xref="paper",
                    yref="paper"
                )
            
        # Update layout for multiple subplots
        fig.update_layout(
            title="Protein Domains",
            xaxis_title='Residual',
            yaxis_title='Domains',
            height=5600 + 100 * len(self.systems),  # Adjust height based on number of systems
            width=1800,  # Adjust width as needed
            autosize=True,
            margin=dict(l=50, r=50, t=50, b=50),  # Adjust margins as needed
            showlegend=False,
            **{f'xaxis{i+1}': dict(
                showline=True,
                showgrid=True,
                showticklabels=True,
                zeroline=False
            ) for i in range(len(self.systems))}
        )

        # Save to HTML
        pio.write_html(fig, file='plots/plot2.html', auto_open=True)

class System(Family):
    def __init__(self, fam_id, sys_id, sys_len, domains):
        self.fam_id = fam_id
        self.sys_id = sys_id
        self.sys_len = sys_len
        self.domains = domains




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

        # Run the CUT command on the output of the CLEAN command
        # cut_process = subprocess.run(CUT_COMMAND, input=clean_process.stdout, capture_output=True, text=True)

        # if cut_process.returncode != 0:
        #     raise Exception(f"Error in CUT command: {cut_process.stderr}")

        # Convert output to DataFrame
        from io import StringIO
        cleaned_data = StringIO(clean_process.stdout)
        df = pd.read_csv(cleaned_data, sep='\t', header=None)

        # Set column names
        header = header.split(", ")
        # header = header[0:5] + header[6:8]
        df.columns = header

        return df  # Return the final output
    
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
                res.append([start, i, -1])  # Use i+1 to account for 1-based index
                start = None

    if start is not None:
        res.append([start, p_length, -1])
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

import matplotlib.colors as mcolors

def get_palette(data):
    uniq_doms = data[DOM_ID].unique()
    # Generate the color palette
    palette = sns.color_palette("husl", n_colors=len(uniq_doms))
    # Convert each color to RGB
    rgb_palette = {domain: mcolors.to_hex(palette[i]) for i, domain in enumerate(uniq_doms)}
    return rgb_palette