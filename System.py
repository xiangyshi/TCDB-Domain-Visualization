'''
Module Name: System.py

Description:
    This module contains the System class that represents a single protein system
    within a TCDB family. It manages the protein's domains and inter-domain regions,
    handling their organization, identification, and analysis.

Dependencies:
    - pandas: For data manipulation
    - numpy: For numerical operations
    - matplotlib: For plotting
    - Domain: For domain representation
    - Hole: For inter-domain region representation
    - util: For utility functions
'''

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
from Family import *
from Hole import *
from Domain import *
import utility.util as util

class System:
    """
    A class representing a single protein system and its domain architecture.

    This class manages a protein's domains and inter-domain regions ("holes"),
    providing methods for their identification and analysis.

    Attributes:
        fam_id (str): TCDB family identifier
        sys_id (str): Unique system identifier
        sys_len (int): Length of the protein sequence
        domains (list[Domain]): List of Domain objects in order of appearance
        accession (str): Protein accession number
        sequence (str): Amino acid sequence
        holes (list[Hole]): List of inter-domain regions
    """

    def __init__(self, fam_id, sys_id, sys_len, domains, sequence, accession):
        """
        Initialize a System object.

        Args:
            fam_id (str): TCDB family identifier
            sys_id (str): Unique system identifier
            sys_len (int): Length of the protein sequence
            domains (list): List of domain tuples (start, end, dom_id, bitscore/evalue)
            sequence (str): Amino acid sequence
            accession (str): Protein accession number
        """
        self.fam_id = fam_id
        self.sys_id = sys_id
        self.sys_len = sys_len
        self.domains = []
        for dom in domains:
            # Check if is rescued domain (only evalue is available)
            if type(dom[3]) == float:
                self.domains.append(Domain(dom[2], dom[0], dom[1], -1, dom[2], (dom[1] - dom[0]) / sys_len, evalue=dom[3]))
            else:
                self.domains.append(Domain(dom[2], dom[0], dom[1], dom[3], dom[2], (dom[1] - dom[0]) / sys_len))
        self.domain_map = {}
        for dom in self.domains:
            if dom.type != "hole":
                if dom.dom_id not in self.domain_map:
                    self.domain_map[dom.dom_id] = []
                self.domain_map[dom.dom_id].append(dom)
        self.accession = accession
        self.sequence = sequence
        self.get_holes()
    
    def check_char(self, char_domains):
        """
        Identifies characteristic domains in the system.

        Args:
            char_domains (list): List of domain IDs considered characteristic

        Returns:
            bool: True if system contains any characteristic domains
        """
        flag = False
        for dom in self.domains:
            if dom.dom_id in char_domains:
                dom.type = "char"
                flag = True
        return flag

    def fill_holes(self):
        # Filter out existing holes and sort domains by start position
        non_hole_domains = [dom for dom in self.domains if dom.type != "hole"]
        sorted_domains = sorted(non_hole_domains, key=lambda x: x.start)
        
        # Create a new list for all domains including holes
        new_domains = []
        
        # Check for gap at the beginning
        if not sorted_domains or sorted_domains[0].start > 0:
            start = 0
            end = sorted_domains[0].start - 1 if sorted_domains else self.sys_len - 1
            hole = Domain("-1", start, end, -1, "-1", (end - start) / self.sys_len)
            new_domains.append(hole)
        
        # Add first domain
        if sorted_domains:
            new_domains.append(sorted_domains[0])
        
        # Check for gaps between domains
        for i in range(1, len(sorted_domains)):
            prev_end = sorted_domains[i-1].end
            curr_start = sorted_domains[i].start
            
            # If there's a gap between the domains
            if prev_end + 1 < curr_start:
                hole = Domain("-1", prev_end + 1, curr_start - 1, -1, "-1", 
                             (curr_start - prev_end - 1) / self.sys_len)
                new_domains.append(hole)
            
            # Add current domain
            new_domains.append(sorted_domains[i])
        
        # Check for gap at the end
        if sorted_domains and sorted_domains[-1].end < self.sys_len - 1:
            start = sorted_domains[-1].end + 1
            end = self.sys_len - 1
            hole = Domain("-1", start, end, -1, "-1", (end - start) / self.sys_len)
            new_domains.append(hole)
        
        # Update domains
        self.domains = new_domains

    def get_holes(self, thresh=50, margin=10):
        """
        Identifies inter-domain regions ("holes") in the system.

        Analyzes the protein sequence to find regions between domains that are
        longer than the threshold length. For each hole, identifies the flanking
        domains within the specified margin.

        Args:
            thresh (int, optional): Minimum length for a region to be considered a hole. 
                                  Defaults to 50.
            margin (int, optional): Number of positions to look for flanking domains. 
                                  Defaults to 10.
        """
        self.holes = []
        for i, dom in enumerate(self.domains):
            if dom.type == "hole" and dom.end - dom.start >= thresh:
                # Find domains within margin distance of hole boundaries
                left_doms, right_doms = util.find_margins(self.domains, 
                                                        dom.start - margin, 
                                                        dom.end + margin)
                ref_doms = []
                names = set()

                # Handle cases based on presence of flanking domains
                if len(left_doms) == 0:
                    for rdom in right_doms:
                        names.add((None, rdom.dom_id))
                        ref_doms.append((None, rdom))
                elif len(right_doms) == 0:
                    for ldom in left_doms:
                        names.add((ldom.dom_id, None))
                        ref_doms.append((ldom, None))
                else:
                    for ldom in left_doms:
                        for rdom in right_doms:
                            names.add((ldom.dom_id, rdom.dom_id))
                            ref_doms.append((ldom, rdom))

                # Create new Hole object
                self.holes.append(Hole(self.sys_id.split('-')[0], i, names, 
                                     ref_doms, dom.start, dom.end, 
                                     self.sequence[dom.start - 1: dom.end]))
    
    def domain_to_list(self):
        """
        Converts domain information to list format, excluding holes.

        Returns:
            list: List of domain tuples (dom_id, start, end, bitscore)
        """
        return [dom.to_tuple() for dom in self.domains if dom.type != "hole"]

    def holes_to_list(self):
        """
        Converts hole information to list format.

        Returns:
            list: List of hole tuples (name, start, end)
        """
        return [hole.to_tuple() for hole in self.holes]

    def generate_csv_row(self):
        """
        Generates a CSV row containing system information.

        Returns:
            list: [accession, length, family_id, subfamily, domains_list, holes_list]
        """
        subfamily = ".".join(self.sys_id.split(".")[:4])
        return [self.accession, self.sys_len, self.fam_id, subfamily, 
                str(self.domain_to_list()), str(self.holes_to_list())]