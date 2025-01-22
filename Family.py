'''
Module Name: Family.py

Description:
    This module contains the Family class that simulates the structure of a family.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from collections import Counter
from config import *
import scipy.stats as stats
from System import *
import util

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
        self.sequence_map = self.get_sequences()
        self.systems = self.get_systems()
        self.char_domains = self.get_char_domains()
        bad_sys = []
        for sys in self.systems:
            sys.check_char(self.char_domains)
            if not sys.has_char:
                bad_sys.append(sys.sys_id)
        #print("Creating palette.")

        self.holes = self.generate_checklist()
        self.gen_palette, self.char_palette, self.hole_palette = self.get_palettes()
    
    def get_sequences(self):
        sequence_map = {}
        faa_path = "./sequences/tcdb-" + self.fam_id + ".faa"
        with open(faa_path, 'r') as sfile:
            lines = sfile.readlines()
            for i in range(len(lines)):
                if i % 2 == 0:
                    sequence_map[lines[i][1:].strip()] = lines[i+1].strip()
        return sequence_map

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
                domains = util.merge_domains(domains)
            res.append(System(self.fam_id, sys_id, sys_len, domains, self.sequence_map[sys_id]))
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
        links = util.find_holes(sys_len, domains_loc)

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

        hole_palette = sns.color_palette("husl", n_colors=len(self.holes))
        hole_palette = {i: mcolors.to_hex(hole_palette[i]) for i in range(len(self.holes))}

        return gen_palette, char_palette, hole_palette

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
        util.combine_svgs(svgs, "general/"  + self.fam_id + ".html")
    
    def plot_char(self):
        svgs = []
        for i, sys in enumerate(self.systems):
            char_doms = [dom for dom in sys.domains if dom[-1] in self.char_domains]
            fig, ax = plt.subplots(figsize=(16, 0.3 * (len(char_doms) + 2)))  # Adjust size as needed
            
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
        util.combine_svgs(svgs, '/char/' + self.fam_id + "-char.html")

    def plot_summary(self):
        fig, ax = plt.subplots(figsize=(16, 0.3 * (len(self.systems) + 2)))
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
            start = util.confidence_interval_mean([i[0] for i in intervals])
            end = util.confidence_interval_mean([i[1] for i in intervals])
            
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
        fig, ax = plt.subplots(figsize=(16, 0.25 * (len(self.char_domains) + 15)))
        for char_dom, intervals in char_domain_map.items():
            for interval in intervals:
                space = np.linspace(interval[0], interval[1], 2)
                ax.plot(space, [cnt] * 2, color=self.char_palette[char_dom.split(' ')[0]], label=char_dom.split(' ')[0], linewidth=8)
            cnt += 1
        
        ax.set_yticks(range(len(char_domain_map)))  # Ensure y-ticks match the number of systems
        ax.set_yticklabels(list(char_domain_map.keys()))  # Set the custom y-tick labels
        
        ax.set_xlabel('Residual Percentile %')
        ax.set_title(self.fam_id + ' Architecture')

        ax.set_xlim(-1, 101)

        svg = "plots/arch/"  + self.fam_id + "-arch.svg"
        fig.savefig(svg)
        plt.close()
        util.combine_svgs([svg], '/arch/' + self.fam_id + "-arch.html")
    
    def plot_holes(self):
        fig, ax = plt.subplots(figsize=(16, 0.3 * (len(self.systems) + 2)))
    
        # Set y-ticks and y-tick labels
        ax.set_yticks(range(len(self.systems)))  # Ensure y-ticks match the number of systems
        sys_ids = {sys.sys_id.split('-')[0]:i for i, sys in enumerate(self.systems)}  # Collect the system IDs
        ax.set_yticklabels(list(sys_ids.keys()))  # Set the custom y-tick labels
        
        for sys in self.systems:
            ax.barh(sys_ids[sys.sys_id.split("-")[0]], sys.sys_len, left=0, color="red", height=0.03)
        
        for i, group in enumerate(self.holes):
            for hole in group:
                ax.barh(sys_ids[hole.sys_id], hole.end - hole.start, color=self.hole_palette[i], left=hole.start, height=0.4)

        # Set axis labels and title
        ax.set_xlabel('Residual')
        ax.set_title(self.fam_id + ' Holes')

        # Save the plot to an HTML file
        fig.savefig("plots/holes/"  + self.fam_id + "-holes.svg")
        plt.close()
    
    def generate_checklist(self):
        holes = []
        pairs = []
        for sys in self.systems:
            for hole in sys.holes:
                holes.append(hole)
        for i in range(len(holes)):
            for j in range(i + 1, len(holes)):
                if util.compare_reference(holes[i], holes[j]):
                    pairs.append([i, j])
        idx_sets = util.get_connected_components(len(holes), pairs)
        
        res = [[holes[i] for i in idx_set] for idx_set in idx_sets]
        for i, idx_set in enumerate(res):
            if len(idx_set) == 1:
                continue
            with open("./holes/" + self.fam_id + "_" + str(i) +  "_holes.fasta", "w") as file:
                for j, hole in enumerate(idx_set):
                    file.write(f">{hole.sys_id}_g{i}_{j}\n")
                    file.write(hole.sequence)
                    file.write("\n")
            
        return res