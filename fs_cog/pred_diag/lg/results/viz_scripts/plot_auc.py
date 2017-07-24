import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from custom import utils

#sns.set_style("whitegrid")
sns.set(font_scale=2, style="whitegrid")

# violin plot for all of the models in the respath variable
respath = {
"lg_XBA.pkl":0,
"lg_XBC.pkl":0,
"lg_XCC.pkl":0,
"lg_XBCCR.pkl":0,
"lg_XBCR.pkl":0,
"lg_XCCR.pkl":0,
"lg_XBCOV.pkl":0,
"lg_XB.pkl":0,
"lg_XC.pkl":0
}
path = "/storage/gablab001/data/genus/fs_cog/pred_diag/lg/results"

for key, val in respath.items():
    respath[key] = os.path.join(path, key)

auc_df = pd.DataFrame(columns=[key.split('_')[1].replace('.pkl', '') for key in respath.keys()])

for key, val in respath.items():
    col = key.split('_')[1].replace('.pkl', '')
    data = utils.read_pickle(val)['auc']
    auc_df[col] = [i[1] for i in data]

order = ['XBA','XBC','XBCCR','XBCOV','XBCR','XB','XCC','XCCR','XC']
sns.violinplot(auc_df, order=order)
plt.title("Distributions of AUC")
fig = plt.gcf()
fig.set_size_inches(12, 7)
fig.savefig('lg_auc_dist.png', dpi=400, bbox_inches="tight")
plt.close()
