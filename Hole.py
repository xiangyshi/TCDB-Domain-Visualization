'''
Module Name: Hole.py

Description:
    This module contains the Hole class that simulates the structure of a hole.
'''
class Hole:
    def __init__(self, sys_id, pos, names, start, end, sequence):
        self.sys_id = sys_id
        self.pos = pos
        self.names = names
        self.start = start
        self.end = end
        self.sequence = sequence
