import os
import numpy as np
import pandas as pd
from custom import utils
from collections import Counter

data_path = "/storage/gablab001/data/genus/fs_cog"
data_name = "FS_COGSCORES.csv"
data_save_path = "/storage/gablab001/data/genus/fs_cog/pred_diag/data_sets"

patient = ['2','Bipolar Disorder','Major Depressive Disorder','Schziophrenia']
control = ['1', 'Control']
cog_scores_to_use = ['SOP_DOMAIN', 'VLM_DOMAIN', 'RPS_DOMAIN']
meta = ['SEX_NEURO', 'AGE_NEURO', 'HANDED', 'STUDY', 'GROUP']
bcvar = ["EstimatedTotalIntraCranialVol"]
cols = colheads.tolist() + cog_scores_to_use + meta + bcvar

idx = lambda t, d: reduce(np.logical_or,
                   [[i == x for x in d.GROUP.values] for i in t])

def load_n_idx(path, name):
    df = pd.read_csv(os.path.join(path, name), low_memory=False)
    pidx, cidx = (idx(patient, df), idx(control, df))
    controls = df.iloc[cidx, :].reset_index(drop=True)
    patients = df.iloc[pidx, :].reset_index(drop=True)
    return pd.concat([controls, patients], axis=0).drop_duplicates('IID')

text_dir = '/storage/gablab001/data/genus/current/genus/text_files_for_indexing'
brain_dir = 'storage/gablab001/data/genus/current/structured/brain/'
thickness = np.genfromtxt(os.path.join(text_dir, '170_columns.txt'), dtype=str)
colheads = thickness

XDF = load_n_idx(data_path, data_name)[cols].dropna().reset_index(drop=True)
y = (~np.logical_or(XDF.GROUP == '1', XDF.GROUP == 'Control')).astype(int).values
assert np.any(~pd.isnull(XDF)), "null values in the data"
site = utils.make_non_singular(utils.encoder(XDF.STUDY.values))
site = pd.DataFrame(data=site,
    columns=['site{}'.format(i) for i in range(site.shape[1])])
handed = utils.make_non_singular(utils.encoder(XDF.HANDED.values))
handed = pd.DataFrame(data=handed,
    columns=['handed{}'.format(i) for i in range(handed.shape[1])])
sex = XDF.SEX_NEURO
age = XDF.AGE_NEURO
icv = XDF[bcvar]
cvars = pd.concat([site, handed, sex, age, icv], axis=1)

XBC = XDF[colheads.tolist() + cog_scores_to_use]
XBA = pd.concat([XBC, cvars], axis=1)
XB = XDF[colheads.tolist()]
XC = XDF[cog_scores_to_use]
XBCOV = pd.concat([XB, cvars], axis=1)
XCC = pd.concat([XC, cvars], axis=1)
XBCR = pd.DataFrame(utils.proj(XB.values, cvars.values), columns=XBC.columns)
XCCR = pd.DataFrame(utils.proj(XC.values, cvars.values), columns=XC.columns)
XBCCR = pd.DataFrame(utils.proj(XBC.values, cvars.values), columns=XBCCR.columns)

for key, val in {'XBC': XBC, 'XBA': XBA, 'XB': XB, 'XC': XC,
                  'XBCOV': XBCOV, 'XCC': XCC, 'XBCR': XBCR,
                  'XCCR': XCCR, 'XBCCR': XBCCR}.items():
    save_name = key + ".csv"
    save_path = os.path.join(data_save_path, save_name)
    val.to_csv(save_path, index=None)

res_save = os.path.join(data_save_path, "response.txt")
np.savetxt(res_save, y, fmt="%i")
