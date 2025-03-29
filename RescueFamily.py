import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from config import *
from Family import Family
import util

class RescueFamily(Family):

    def __init__(self, rescue_data, fam_id):
        self.domain_priority = []
        super().__init__(rescue_data, fam_id)
    
    def plot_char_rescue(self):
            svgs = []
            colors = ["green", "orange", "red"]
            for i, sys in enumerate(self.systems):
                char_doms = [dom for dom in sys.domains if dom.dom_id in self.char_domains]
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
                    curr_sys = self.data[self.data[SYS_ID] == sys.sys_id]
                    rround = curr_sys[curr_sys[DOM_ID] == dom.dom_id]["rescue round"].iloc[0]
                    space = np.linspace(dom.start, dom.end, 2)
                    ax.plot(space, [cnt] * 2, color=colors[int(rround)], label=dom.dom_id, linewidth=8)

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
            util.combine_svgs(svgs, '/resc/' + self.fam_id + "-resc.html")