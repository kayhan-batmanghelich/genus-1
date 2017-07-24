import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from custom import utils
from collections import Counter
import matplotlib.pyplot as plt

# command line inputs
result_name = sys.argv[1]
header_name = sys.argv[2]
save_name = sys.argv[3]

# headers
hp = "/storage/gablab001/data/genus/fs_cog/pred_diag/column_headers"

with open(os.path.join(hp, header_name), "r") as tmp:
    headers = np.array(tmp.read().strip("\n").split("\n")) 

coefs = utils.read_pickle("../{}".format(result_name))['coef']
nz = lambda x: np.nonzero(x)[0]
cog_scores = ['SOP','RPS','VLM']
cols = []

for coef in coefs:
    selected = headers[nz(coef[1])]
    cols.extend(selected.tolist())

counts = dict(Counter(cols))
to_keep = {}

for key, val in counts.items():
    s = 'DOMAIN'
    if ('SOP_'+s in key) or ('VLM_'+s in key) or ('RPS_'+s in key):
        to_keep[key] = val

df = pd.DataFrame({'Cog Variable': [key for key in to_keep.keys()], 
                   'Times Selected': [val for val in to_keep.values()]})

sns.set(style="whitegrid", font_scale=2) 
sns.barplot(x='Cog Variable', y='Times Selected', data=df)
plt.ylabel("Times Selected")
save_name = "../images/{}_cog_scores.png".format(save_name)
plt.savefig(save_name, bbox_inches="tight", dpi=300)
plt.close()
   
