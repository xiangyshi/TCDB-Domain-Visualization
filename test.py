from util import *
import seaborn as sns

input = 'outFile.cdd'
df = get_clean(input)

uniq_doms = df[DOM_ID].unique()
palette = sns.color_palette("husl", n_colors=len(uniq_doms))
domain_to_color = {domain: palette[i] for i, domain in enumerate(uniq_doms)}