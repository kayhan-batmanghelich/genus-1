import os
import numpy as np
import pandas as pd
from custom import utils
"""
cvars.csv  response.txt  XBA.csv  XBCCR.csv  XBC.csv  
XBCOV.csv  XBCR.csv  XB.csv  XCC.csv  XCCR.csv  XC.csv  XDF.csv
"""
datasets_path = "/storage/gablab001/data/genus/fs_cog/exp1/singularity/datasets_for_r"
#cvars = pd.read_csv(os.path.join(datasets_path, "cvars.csv"))
XBA = pd.read_csv(os.path.join(datasets_path, "XBA.csv"))
site_cols = np.array(['site{}'.format(i) for i in range(17)])
sites = (XBA[site_cols]*np.arange(17)).sum(1).values
XBA['sites_combined'] = sites
minus_ohe_sites = np.setdiff1d(XBA.columns.values, site_cols)
XBA = XBA[minus_ohe_sites]
data_by_site = XBA.groupby('sites_combined')


from sklearn import linear_model
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import StratifiedKFold
