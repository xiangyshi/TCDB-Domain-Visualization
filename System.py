import subprocess
import platform
import os
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
import util

class System:
    def __init__(self, fam_id, sys_id, sys_len, domains, sequence):
        self.fam_id = fam_id
        self.sys_id = sys_id
        self.sys_len = sys_len
        self.domains = domains
        self.sequence = sequence
        self.get_holes()
    
    def check_char(self, char_domains):
        for dom in self.domains:
            if dom[2] in char_domains:
                self.has_char = True
                return
        self.has_char = False
    
    def get_holes(self, thresh=50, margin=10):
        self.holes = []
        for i, dom in enumerate(self.domains):
            if dom[-1] == "-1" and dom[1] - dom[0] >= thresh:
                left_doms, right_doms = util.find_margins(self.domains, dom[0] - margin, dom[1] + margin)

                names = set()
                if len(left_doms) == 0:
                    for rdom in right_doms:
                        names.add((None, rdom[-1]))
                elif len(right_doms) == 0:
                    for ldom in left_doms:
                        names.add((ldom[-1], None))
                else:
                    for ldom in left_doms:
                        for rdom in right_doms:
                            names.add((ldom[-1], rdom[-1]))

                self.holes.append(Hole(self.sys_id.split('-')[0], i, names, dom[0], dom[1], self.sequence[dom[0] - 1: dom[1]]))
