import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from config import *
from Family import Family
import util
import networkx as nx

class RescueFamily(Family):

    def __init__(self, rescue_data, fam_id):
        super().__init__(rescue_data, fam_id)
        self.filter_domains()
    
    def filter_domains(self):
        for sys in self.systems:
            # Create a graph
            G = nx.Graph()

            sys_domains = [dom for dom in sys.domains if dom.type != "hole"]
            # Add nodes for each domain
            for i, dom in enumerate(sys_domains):
                G.add_node(i, domain=dom)

            # Add edges for overlapping domains
            for i in range(len(sys_domains)):
                for j in range(i + 1, len(sys_domains)):
                    if util.is_overlap(sys_domains[i], sys_domains[j], self.data):
                        G.add_edge(i, j)

            # Find connected components
            connected_components = list(nx.connected_components(G))
            print("connected components: ", sys.sys_id, connected_components)
            # Store the overlapping domains for the system
            sys.overlapping_domains = [[G.nodes[i]['domain'] for i in component] for component in connected_components]

            # Select representative domain for each connected component
            sys_new_domains = []
            for component in connected_components:
                # Select the domain with the highest score
                max_score = float('-inf')
                representative_domain = None
                for i in component:
                    
                    domain = G.nodes[i]['domain']
                    score = util.score_domain(domain, sys.sys_len)  # Calculate the score for the domain
                    if score > max_score:
                        max_score = score
                        representative_domain = domain
                
                # Store the representative domain for the component
                if representative_domain:
                    print("representative domain: ", component, representative_domain)
                    sys_new_domains.append(representative_domain)
            sys.domains = sys_new_domains

    def plot_char_rescue(self):
            svgs = []
            colors = ["green", "orange", "red"]
            for i, sys in enumerate(self.systems):
                sys_doms = sys.domains
                fig, ax = plt.subplots(figsize=(16, 0.3 * (len(sys_doms) + 2)))  # Adjust size as needed
                
                ax.set_title(sys.sys_id)
                ax.set_xlabel('Residual')
                ax.set_ylabel('Domains')

                # Plot domains
                space = np.linspace(0, sys.sys_len - 1, 2)
                cnt = 0
                dom_ids = [sys.sys_id.split('-')[0]]
                ax.plot(space, [cnt] * 2)
                
                for dom in sys_doms:
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