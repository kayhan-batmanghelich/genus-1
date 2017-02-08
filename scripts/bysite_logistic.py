from __future__ import division
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn import linear_model
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

def covar_idx(data, covars):
    return [data.columns.get_loc(col) for col in covars]

def encode(data):
    encoded_label = data.name
    le = LabelEncoder()
    ohe = OneHotEncoder()
    temp_encode = le.fit_transform(data)
    final_encode = ohe.fit_transform(temp_encode.reshape(-1,1)).toarray()
    return pd.DataFrame(final_encode, columns = le.classes_)

covars_1 = [
    'SEX', 'HEAD_COIL_NCHANNELS',
    'FS_VERSION', 'MAGN_FIELD_STR',
    'HANDED', 'STUDY', 'IID', 'GROUP',
    'AGE_MRI'
]

def res(y):
    p = 'Schizophrenia'
    return np.array([1. if val == p else 0 for val in y])

data = pd.read_csv('GENUS_FS_ATLAS_D.csv', low_memory = False)
covariates_non_encoded = data.ix[:, covar_idx(data, covars_1)]
data_brain_vals_only = data.ix[:, :-25]
data_to_dropna = pd.concat([data_brain_vals_only, covariates_non_encoded], axis = 1)
dna = data_to_dropna.dropna(axis=0)
dna.index = [x for x in range(dna.shape[0])]
diagnosis_grouping = dna.groupby('GROUP')
scz = diagnosis_grouping.get_group('Schizophrenia')
cnt = diagnosis_grouping.get_group('Control')
diag_scz_cnt = pd.concat([scz, cnt], axis=0)
diag_scz_cnt.index = [x for x in range(diag_scz_cnt.shape[0])]
groups = diag_scz_cnt.groupby('STUDY')

def rem_redudant(encoded):
    to_drop = []
    for col in range(1, encoded.shape[1]):
