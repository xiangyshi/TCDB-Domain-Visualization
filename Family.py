'''
Module Name: Family.py

Description:
    This module contains the Family class that simulates the structure of a protein family.
    It handles the organization of TCDB (Transporter Classification Database) data and 
    provides methods for analyzing and visualizing domain architectures across multiple
    protein systems within a family.

Dependencies:
    - numpy: For numerical operations
    - matplotlib: For plotting and visualization
    - seaborn: For color palettes
    - scipy.stats: For statistical analysis
    - System: For handling individual protein systems
    - util: For utility functions
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from collections import Counter
from utility.config import *
import scipy.stats as stats
from System import *
import utility.util as util

class Family:
    '''
    A class representing a protein family and its domain architecture.

    The Family class organizes protein systems that belong to the same TCDB family,
    analyzes their domain composition, and provides various visualization methods
    for domain architecture analysis.

    Attributes:
        data (pd.DataFrame): Cleaned rpblast output data for this family
        fam_id (str): Unique TCDB accession of this family
        merge (bool): Whether to merge overlapping domains
        sequence_map (dict): Mapping of system IDs to their sequences
        systems (list): Collection of System objects in the family
        char_domains (list): Set of characteristic domains in this family
        holes (list): List of inter-domain regions
        gen_palette (dict): Color palette for general domains
        char_palette (dict): Color palette for characteristic domains
        hole_palette (dict): Color palette for holes
        max_sys_len (int): Length of the longest protein in the family
    '''
     
    def __init__(self, data, fam_id, output, merge=True):
        '''
        Initialize a Family object.

        Args:
            data (pd.DataFrame): Contains the cleaned output of this family from rpblast
            fam_id (str): Unique TCDB accession of this family
            merge (bool, optional): Whether to merge overlapping domains. Defaults to True
        '''
        self.data = data
        self.fam_id = fam_id
        self.merge = merge
        self.sequence_map = self.get_sequences()
        self.systems = self.get_systems()
        self.char_domains = self.get_char_domains()
        bad_sys = []
        for sys in self.systems:
            has_char = sys.check_char(self.char_domains)
            if not has_char:
                bad_sys.append(sys.sys_id)
                
        self.holes = self.generate_checklist()
        self.gen_palette, self.char_palette, self.hole_palette = self.get_palettes()
        self.output = output

        # Create and clean output directory
        os.makedirs(os.path.join(self.output, self.fam_id), exist_ok=True)
        
    def get_sequences(self):
        '''
        Reads and maps protein sequences from FASTA file.

        Returns:
            dict: Mapping of system IDs to their protein sequences
        '''
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
        Constructs System objects for each unique system in the family data.
        
        Returns:
            list[System]: List of constructed System objects
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
            res.append(System(self.fam_id, sys_id, sys_len, domains, self.sequence_map[sys_id], sys_id.split('-')[-1]))
        return res
    
    def get_domains(self, sys_len, data):
        '''
        Retrieves domain information for a specific system.
        
        Args:
            sys_len (int): Total length of the system protein
            data (pd.DataFrame): Data containing hits only within this system

        Returns:
            list: Domain hits and inter-domain regions, sorted by position
        '''
        try:
            domains = data[[DOM_START, DOM_END, DOM_ID, BIT_SCORE]].sort_values(DOM_START)
        except:
            domains = data[[DOM_START, DOM_END, DOM_ID, E_VALUE]].sort_values(DOM_START)
        domains_loc = np.array(domains[[DOM_START, DOM_END]]).tolist()
        links = util.find_holes(sys_len, domains_loc)

        domains = links + np.array(domains).tolist()
        domains = sorted(domains, key=lambda x: x[0])
        return domains
    
    def get_char_domains(self, thresh=0.5):
        '''
        Identifies characteristic domains of the family.

        A domain is considered characteristic if it appears in more than
        the threshold percentage of systems in the family.

        Args:
            thresh (float): Threshold presence percentage (0-1) for characteristic domains

        Returns:
            list: Domain IDs considered characteristic for this family
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
        Creates color palettes for domain visualization.

        Returns:
            tuple: Contains (general_domain_palette, characteristic_domain_palette, hole_palette)
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
        Generates domain architecture plots for each system.

        Args:
            mode (str): Plot mode - 'char' prioritizes characteristic domains in visualization
        '''
        svgs = []
        for i, sys in enumerate(self.systems):
            size = len([0 for dom in sys.domains if dom.type != "hole"])
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
                gen_doms = [dom for dom in sys.domains if dom.type == "char"] + [dom for dom in sys.domains if dom.type != "char"]

            for dom in gen_doms:
                if dom.type == "hole":
                    continue
                cnt -= 1
                dom_ids.append(dom.dom_id)
                space = np.linspace(dom.start, dom.end, 2)
                if dom.dom_id in self.gen_palette:
                    ax.plot(space, [cnt] * 2, color=self.gen_palette[dom.dom_id], label=dom.dom_id, linewidth=8)
                else:
                    ax.plot(space, [cnt] * 2, color=CHAR_COLOR, label=dom.dom_id, linewidth=8)

            # Set consistent axis limits
            ax.set_xlim(0, self.max_sys_len)  # Adjust based on your data
            ax.set_ylim(cnt - 1, 1)
            
            ax.set_yticks(range(0, cnt - 1, -1))
            ax.set_yticklabels(dom_ids)

            # Save as SVG
            svg_file = os.path.join(self.output, f"{self.fam_id}/general-{sys.sys_id}.svg")
            fig.savefig(svg_file, format='svg', bbox_inches='tight', pad_inches=0.2)
            svgs.append(svg_file)
            plt.close(fig)  # Close the figure to avoid memory issues
        util.combine_svgs(svgs, f"general-{self.fam_id}.html", os.path.join(self.output, self.fam_id))
    
    def plot_char(self):
        '''
        Generates plots showing only characteristic domains for each system.
        '''
        svgs = []
        for i, sys in enumerate(self.systems):
            char_doms = [dom for dom in sys.domains if dom.type == "char"]
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
                dom_ids.append(dom.dom_id)
                space = np.linspace(dom.start, dom.end, 2)
                ax.plot(space, [cnt] * 2, color=self.char_palette[dom.dom_id], label=dom.dom_id, linewidth=8)

            # Set consistent axis limits
            ax.set_xlim(0, self.max_sys_len)  # Adjust based on your data
            ax.set_ylim(cnt - 1, 1)
            
            ax.set_yticks(range(0, cnt - 1, -1))
            ax.set_yticklabels(dom_ids)

            # Save as SVG
            svg_file = os.path.join(self.output, f"{self.fam_id}/char-{sys.sys_id}.svg")
            fig.savefig(svg_file, format='svg', bbox_inches='tight', pad_inches=0.2)
            svgs.append(svg_file)

            plt.close(fig)  # Close the figure to avoid memory issues
        util.combine_svgs(svgs, f"char-{self.fam_id}.html", os.path.join(self.output, self.fam_id))

    def plot_summary(self):
        '''
        Generates a summary plot showing domain coverage across all systems.
        '''
        fig, ax = plt.subplots(figsize=(16, 0.3 * (len(self.systems) + 2)))
        for i, sys in enumerate(self.systems):
            for dom in sys.domains:
                if dom.type == "hole":
                    ax.barh(i, dom.end - dom.start, left=dom.start, height=0.03, color='blue')
                else:
                    ax.barh(i, dom.end - dom.start, left=dom.start, height=0.4, color='blue')

        # Set y-ticks and y-tick labels
        ax.set_yticks(range(len(self.systems)))  # Ensure y-ticks match the number of systems
        sys_ids = [sys.sys_id.split('-')[0] for sys in self.systems]  # Collect the system IDs
        ax.set_yticklabels(sys_ids)  # Set the custom y-tick labels
        
        # Set axis labels and title
        ax.set_xlabel('Residual')
        ax.set_title(self.fam_id + ' Summary')

        # Save the plot to an SVG file
        svg_file = os.path.join(self.output, f"{self.fam_id}/summary-{self.fam_id}.svg")
        fig.savefig(svg_file)
        plt.close()
        util.combine_svgs([svg_file], f"/summary-{self.fam_id}.html", os.path.join(self.output, self.fam_id))

    def plot_arch(self):
        '''
        Generates a consensus architecture plot for the family.
        
        Analyzes the positions of characteristic domains across all systems
        to create a representative architecture for the family. Requires at
        least 3 systems for meaningful analysis.
        '''
        if len(self.systems) <= 2:
            print("System count too small, unable to accurately construct architecture. Skipping arch plot...")
            return
        
        char_domain_locs = {k:[] for k in self.char_domains}
        char_domain_map = {}
        # Store char domain locations in format "char_dom_id" : [[#char_dom_id appearances in sys1], [#char_dom_id appearances in sys2], ... , [#char_dom_id appearances in sysn]]
        for system in self.systems:
            sys_dom_locs = {k:[] for k in self.char_domains}
            for dom in system.domains:
                if dom.dom_id in char_domain_locs:
                    sys_dom_locs[dom.dom_id].append([dom.start / system.sys_len, dom.end / system.sys_len])
            
            
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

        svg = os.path.join(self.output, f"{self.fam_id}/arch-{self.fam_id}.svg")
        fig.savefig(svg)
        plt.close()
        util.combine_svgs([svg], f"/arch-{self.fam_id}.html", os.path.join(self.output, self.fam_id))
    
    def plot_holes(self):
        '''
        Generates visualization of inter-domain regions (holes) across systems.
        '''
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
        svg = os.path.join(self.output, f"{self.fam_id}/holes-{self.fam_id}.svg")
        fig.savefig(svg)
        plt.close()
        util.combine_svgs([svg], f"/holes-{self.fam_id}.html", os.path.join(self.output, self.fam_id))
    
    def generate_checklist(self):
        '''
        Identifies and groups similar inter-domain regions across systems.

        Returns:
            list: Groups of similar holes across systems
        '''
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
    
    def generate_csv_rows(self):
        '''
        Generates CSV output rows for all systems in the family.

        Returns:
            list: CSV rows containing system information
        '''
        rows = []
        for sys in self.systems:
            rows.append(sys.generate_csv_row())
        return rows

    def plot(self, options=["general", "char", "arch", "holes", "summary"]):
        if "general" in options:
            self.plot_general()
        if "char" in options:
            self.plot_char()
        if "arch" in options:
            self.plot_arch()
        if "holes" in options:
            self.plot_holes()
        if "summary" in options:
            self.plot_summary()