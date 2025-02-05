'''
Module Name: Domain.py

Description:
    This module contains the Domain class that represents a protein domain or
    inter-domain region ("hole") within a protein sequence. Each domain has
    a unique identifier, start/end positions, bitscore, and type designation.
'''

class Domain:
    """
    A class representing a protein domain or inter-domain region.

    Attributes:
        dom_id (str): Unique identifier for the domain
        start (int): Start position of the domain in the protein sequence
        end (int): End position of the domain in the protein sequence
        bitscore (float): Bitscore from domain alignment
        type (str): Type of region - either "dom" for domain or "hole" for inter-domain region
    """

    def __init__(self, dom_id, start, end, bitscore, type):
        """
        Initialize a new Domain instance.

        Args:
            dom_id (str): Unique identifier for the domain
            start (int): Start position in protein sequence
            end (int): End position in protein sequence
            bitscore (float): Bitscore from domain alignment
            type (str): "-1" for holes, any other value for domains
        """
        self.dom_id = dom_id
        self.start = start
        self.end = end
        self.bitscore = bitscore
        self.type = "dom" if type != "-1" else "hole"

    def to_tuple(self):
        """
        Convert domain information to a tuple format.

        Returns:
            tuple: Contains (dom_id, start, end, bitscore)
        """
        return (self.dom_id, self.start, self.end, self.bitscore)