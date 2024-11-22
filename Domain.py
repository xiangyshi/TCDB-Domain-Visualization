'''
Module Name: Domain.py

Description:
    This module contains the Domain class that simulates the structure of a domain.
'''

class Domain:
    def __init__(self, dom_id, start, end, type):
        self.dom_id = dom_id
        self.start = start
        self.end = end
        self.type = type