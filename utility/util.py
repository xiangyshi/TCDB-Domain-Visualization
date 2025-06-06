'''
Module Name: util.py

Description:
    This module provides utility functions for domain architecture analysis,
    supporting the model classes (Family, System, and Hole). Functions include
    data cleaning, domain analysis, statistical operations, and visualization helpers.

Dependencies:
    - subprocess: For running system commands
    - platform: For system architecture detection
    - pandas: For data manipulation
    - numpy: For numerical operations
    - scipy.stats: For statistical analysis
    - networkx: For graph operations
    - config: For configuration settings
'''

import os

import numpy as np
from utility.config import *
import scipy.stats as stats
import networkx as nx
from Domain import *


def find_holes(p_length, domain_regions):
    """
    Identifies inter-domain regions ("holes") in a protein sequence.

    Creates a binary mask of the protein sequence where 1 indicates
    positions not covered by any domain, then identifies continuous
    regions of uncovered positions.

    Args:
        p_length (int): Length of the protein sequence
        domain_regions (list): List of domain regions as [start, end] pairs

    Returns:
        list: List of holes as [start, end, "-1", "-1"] where -1 indicates
              hole type and placeholder bitscore
    """
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
                res.append([start, i, "-1", "-1"])  # Use i+1 to account for 1-based index
                start = None

    if start is not None:
        res.append([start, p_length, "-1", "-1"])
    return res

def combine_svgs(svgs, filename, output):
    """
    Combines multiple SVG files into a single HTML file.

    Creates an HTML document containing all provided SVG content
    and removes the original SVG files.

    Args:
        svgs (list): List of paths to SVG files
        filename (str): Output HTML filename

    Returns:
        None
    """
    html_content = ''
    for svg in svgs:
        with open(svg, 'r') as file:
            svg_content = file.read().strip()
            # Ensure the SVG content is well-formed and includes the correct XML namespaces
            html_content += f'<div>{svg_content}</div>'
        os.remove(svg)
    html_content += '</body></html>'
    with open(f'{output}/{filename}', "w") as f:
        f.write(html_content)

def merge_domains(domains):
    """
    Merges overlapping domains that share the same ID.

    Combines adjacent or overlapping domains with the same identifier
    into larger domains, preserving the highest bitscore.

    Args:
        domains (list): List of domain tuples (start, end, id, bitscore)

    Returns:
        list: Sorted list of merged domain tuples
    """

    # Sort intervals by id and then by the start of the interval
    domains = sorted(domains, key=lambda x: (x[2], x[0]))

    merged = []
    current_start, current_end, current_id, current_bitscore = domains[0]

    for i in range(1, len(domains)):
        start, end, id_, bitscore = domains[i]

        if id_ == current_id:
            # If the intervals overlap, merge them
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                # If they don't overlap, add the current interval to the result
                merged.append((current_start, current_end, current_id, current_bitscore))
                current_start, current_end, current_id, current_bitscore = start, end, id_, bitscore
        else:
            # If the ID changes, add the current interval to the result
            merged.append((current_start, current_end, current_id, current_bitscore))
            current_start, current_end, current_id, current_bitscore = start, end, id_, bitscore

    # Add the last interval
    merged.append((current_start, current_end, current_id, current_bitscore))
    return sorted(merged, key=lambda x : x[0])

def confidence_interval_mean(coords, confidence_level=0.95):
    """
    Calculates the confidence interval mean for a set of coordinates.

    For single-value sets, returns the value. For multiple values,
    calculates the confidence interval using Student's t-distribution.

    Args:
        coords (list): List of coordinate values
        confidence_level (float): Confidence level (0 to 1), defaults to 0.95

    Returns:
        float: Mean value of the confidence interval
    """

    if len(set(coords)) == 1:
        return np.mean(coords)
    
    mean = np.mean(coords)
    std_err = stats.sem(coords)
    df = len(coords) - 1

    ci = stats.t.interval(confidence_level, df, loc=mean, scale=std_err)
    return np.mean(ci)

def find_margins(domains, margin_left, margin_right):
    """
    Identifies domains that overlap with specified margin regions.

    Searches for domains that overlap with regions flanking a hole,
    used to identify potential reference domains.

    Args:
        domains (list): List of Domain objects
        margin_left (int): Left margin position
        margin_right (int): Right margin position

    Returns:
        tuple: (left_domains, right_domains) Lists of domains overlapping
               with left and right margins respectively
    """

    left_doms = []
    right_doms = []
    
    for dom in domains:
        if dom.start <= margin_left and margin_left <= dom.end:
            left_doms.append(dom)
        if dom.start <= margin_right and margin_right <= dom.end:
            right_doms.append(dom)

    return left_doms, right_doms

def compare_reference(hole1, hole2):
    """
    Compares two holes to determine if they share reference domains.

    Args:
        hole1 (Hole): First hole object
        hole2 (Hole): Second hole object

    Returns:
        bool: True if holes share at least one reference domain pair
    """
    return len(hole1.names & hole2.names) != 0

def get_connected_components(n, pair_list):
    """
    Groups synonymous holes using a graph-based approach.

    Uses NetworkX to identify connected components in a graph where
    nodes are holes and edges indicate synonymous relationships.

    Args:
        n (int): Total number of holes
        pair_list (list): List of hole pairs that are synonymous

    Returns:
        list: List of sets, each containing indices of synonymous holes
    """
    G = nx.Graph()
    
    # Add all nodes (1-indexed)
    G.add_nodes_from(range(0, n))
    
    # Add all edges (the pairs from the pair_list)
    G.add_edges_from(pair_list)
    
    # Get connected components
    connected_components = list(nx.connected_components(G))
    
    return connected_components

def score_domain(domain, protein_length) -> float:
    # suppress divide by zero warnings
    with np.errstate(divide='ignore'):
        nlog = -np.log(domain.evalue)
        
    # Handle inf values when evalue is 0
    if np.isinf(nlog) or np.isnan(nlog):
        nlog = 1000000
        
    return (domain.end - domain.start / protein_length) * nlog

def is_overlap(dom1, dom2, mutual_overlap=0.2):
    """
    Given two domains and a dataframe containing their information within a system (protein), check if they are overlapping.
    Overlapping is defined by the value of mutual_overlap, which is the overlapping region divided by the longer domain.

    Args:
        dom1 (Domain): The first domain object
        dom2 (Domain): The second domain object
        mutual_overlap (float): Minimum required overlap fraction (0 to 1)

    Returns:
        bool: True if the domains overlap, False otherwise
    """

    # Calculate the overlap
    overlap_start = max(dom1.start, dom2.start)
    overlap_end = min(dom1.end, dom2.end)

    # If there is no overlap
    if overlap_start >= overlap_end:
        return False

    # Calculate the lengths of the domains
    length_dom1 = dom1.end - dom1.start
    length_dom2 = dom2.end - dom2.start

    # Calculate the mutual overlap
    overlap_length = overlap_end - overlap_start
    shorter_length = min(length_dom1, length_dom2)
    # Check if the overlap meets the mutual overlap requirement
    return (overlap_length / shorter_length) >= mutual_overlap
