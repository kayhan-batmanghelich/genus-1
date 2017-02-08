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
        if (encoded.ix[:, col - 1] + encoded.ix[:, col]).mean() == 1:
            to_drop.append(encoded.columns[col])
    return encoded.drop(to_drop,axis=1)

def logistic(X, y, splits):
    out = {'auc':[], 'coef':[]}
    cv_out = StratifiedShuffleSplit(n_splits = splits)
    cv_in = StratifiedShuffleSplit(n_splits = splits - 1)
    for train, test in cv_out.split(X, y):
        clf = Pipeline([
            ('scale', StandardScaler()),
            ('lgr', linear_model.LogisticRegressionCV(
                        cv = cv_in,
                        penalty = 'l1',
                        solver = 'liblinear'))])
        clf.fit(X[train], y[train])
        out['auc'].append(roc_auc_score(y[test], clf.predict(X[test])))
        out['coef'].append(clf.named_steps['lgr'].coef_)
    return np.mean(out['auc']), out['coef']

def reidx(df):
    df1 = df.copy()
    df1.index = [x for x in range(df.shape[0])]
    return df1

results = {}

for name, group in groups:
    temp_encode = []
    cvars = group.ix[:, covar_idx(group, covars_1)]
    for col in range(cvars.shape[1])[:-3]:
        temp_encode.append(encode(cvars.ix[:, col]))
    encoded = pd.concat(temp_encode, axis=1)
    encoded = reidx(rem_redudant(encoded))
    brain_data = reidx(group.ix[:, :-25])
    age = reidx(group[['AGE_MRI']])
    model_data = pd.concat([brain_data, age, encoded], axis=1)
    try:
        results[name] = logistic(X = model_data,
                                 y = ,
                                 splits = 4)
    except ValueError:
        results[name] = 'No variance in population'
