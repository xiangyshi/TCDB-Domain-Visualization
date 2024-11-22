'''
Module Name: util.py

Description:
    This module provides utility functions to the model classes
    such as `Family`, `System`, and `Hole`, as well as cdd extractions
    for the driver code.
'''
import subprocess
import platform
import os
import pandas as pd
import numpy as np
from collections import Counter
from config import *
import scipy.stats as stats
import networkx as nx

def get_clean(in_file) -> pd.DataFrame:
    """
    Cleans the input data by removing comments and selecting specific fields,
    then returns the cleaned data as a pandas DataFrame.

    Args:
        in_file (String): path to input cdd file.

    Returns:
        pd.Dataframe: Dataframe format of the cdd file for further data extraction.
    """

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


def find_holes(p_length, domain_regions):
    """
    Identifies the holes of a protein sequence given the domains

    Args:
        p_length (int): the length of the protein sequence.
        domain_regions (list): a list containing domains (represented in a list with [start, end]).

    Returns:
        list: A list containing the holes of the protein sequence.
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
                res.append([start, i, "-1"])  # Use i+1 to account for 1-based index
                start = None

    if start is not None:
        res.append([start, p_length, "-1"])
    return res

def combine_svgs(svgs, filename):
    """
    Combines a list of svg files into one html file.

    Args:
        svgs (list): contains a list of file paths to svg files.
        filename (String): output destination of html file.

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
    
    with open(f'plots/{filename}', "w") as f:
        f.write(html_content)

def merge_domains(domains):
    """
    Merges overlapping domains (with same id) in to a larger domain.

    Args:
        domains (list): Description of the first parameter, its purpose, and any
            relevant details.

    Returns:
        list: Sorted (based on start index) list of merged domains.
    """

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
    """
    Performs a confidence interval test.

    Args:
        coords (list): A list of coordinate pairs.
        confidence_level (float, 0 to 1): The desired confidence_level.
    Returns:
        list: A pair of coordinates that represents the inferred average
            of the coordinates provided.
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
    Searches for reference domains on the left/right of a hole.
    Function finds any domain that overlapps with the margin segments:
    --------DOMAINS--------|--margin--|-----hole-----|--margin--|--------DOMAINS--------

    Args:
        domains (list): The list of domains of the system containing the hole.
        margin_left (int): start of hole minus margin
        margin_right (int): end of hole plus margin

    Returns:
        tuple: A 2-tuple containing the potential reference domains to the left/right of the hole.
    """

    left_doms = []
    right_doms = []
    
    for dom in domains:
        if dom[0] <= margin_left and margin_left <= dom[1]:
            left_doms.append(dom)
        if dom[0] <= margin_right and margin_right <= dom[1]:
            right_doms.append(dom)

    return left_doms, right_doms

def compare_reference(hole1, hole2):
    """
    Compares two holes to see if they have a common reference domain pair.

    Args:
        hole1 (Hole): First hole to compare.
        hole2 (Hole): Second hole to compare.

    Returns:
        boolean: True if the two holes are synonymous, False otherwise.
    """
    return len(hole1.names & hole2.names) != 0

def get_connected_components(n, pair_list):
    """
    Union-Find / Disjoint-Set algorithm to group synonymous holes.

    Args:
        n (int): the total number of holes in a family.
        pair_list (list): a list containing pairs of synonymous holes.

    Returns:
        list: A list of lists, each sublist containing a group of synonymous hole that
            is not synonmyous with any other group.
    """
    G = nx.Graph()
    
    # Add all nodes (1-indexed)
    G.add_nodes_from(range(0, n))
    
    # Add all edges (the pairs from the pair_list)
    G.add_edges_from(pair_list)
    
    # Get connected components
    connected_components = list(nx.connected_components(G))
    
    return connected_components