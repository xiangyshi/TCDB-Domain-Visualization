'''
Module Name: Hole.py

Description:
    This module contains the Hole class that represents an inter-domain region
    ("hole") between protein domains. A hole is characterized by its position,
    sequence, and the domains that flank it on either side.

Dependencies:
    None
'''

class Hole:
    """
    A class representing an inter-domain region in a protein sequence.

    Attributes:
        sys_id (str): Identifier of the protein system containing this hole
        pos (int): Position of the hole in the sequence of holes
        names (list): List of possible names for this hole
        ref_doms (list): List of domain pairs that flank this hole [(left_dom, right_dom),...]
        start (int): Start position of the hole in the protein sequence
        end (int): End position of the hole in the protein sequence
        sequence (str): Amino acid sequence of the hole region
        best_name (str): Best representative name for the hole based on flanking domains
    """

    def __init__(self, sys_id, pos, names, ref_doms, start, end, sequence):
        """
        Initialize a new Hole instance.

        Args:
            sys_id (str): Identifier of the protein system
            pos (int): Position in the sequence of holes
            names (list): Possible names for this hole
            ref_doms (list): List of domain pairs that flank this hole
            start (int): Start position in protein sequence
            end (int): End position in protein sequence
            sequence (str): Amino acid sequence of the hole
        """
        self.sys_id = sys_id
        self.pos = pos
        self.names = names
        self.ref_doms = ref_doms
        self.start = start
        self.end = end
        self.sequence = sequence
        self.best_name = self.get_best_name()

    def to_tuple(self):
        """
        Convert hole information to a tuple format.

        Returns:
            tuple: Contains (best_name, start, end)
        """
        return (self.best_name, self.start, self.end)
    
    def get_best_name(self):
        """
        Determine the best name for the hole based on its flanking domains.
        
        The name is constructed from the highest scoring domains on either side.
        If no domain exists on either side, "BEGIN" or "END" is used respectively.

        Returns:
            str: Name in format "LEFT_DOMAIN to RIGHT_DOMAIN"
        """
        # Extract left and right domains, filtering out None values
        left_doms = [pair[0] for pair in self.ref_doms if pair[0] is not None]
        right_doms = [pair[1] for pair in self.ref_doms if pair[1] is not None]
        
        # Sort domains by bitscore to get best matches
        left_doms.sort(key=lambda x: x.bitscore, reverse=True)
        right_doms.sort(key=lambda x: x.bitscore, reverse=True)
        
        # Default names if no domains found
        left_best = "BEGIN"
        right_best = "END"
        
        # Get highest scoring domain IDs if available
        try:
            left_best = left_doms[0].dom_id
        except:
            pass
        try:
            right_best = right_doms[0].dom_id
        except:
            pass
        
        return left_best + " to " + right_best
