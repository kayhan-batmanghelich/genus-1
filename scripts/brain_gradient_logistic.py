import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
from custom.custom import encode
from custom.custom import res
from custom.custom import remove_redudant
from custom.custom import reidx

def logistic(X, y, splits):
    results = {'auc':[], 'coef':[]}
    cv_out = StratifiedShuffleSplit(n_splits=splits)
    cv_in = StratifiedShuffleSplit(n_splits=splits - 1)
    for train, test in cv_out.split(X, y):
        clf = Pipeline([
            ('scale', StandardScaler()),
            ('lgr', linear_model.LogisticRegressionCV(
                cv = cv_in,
                penalty = 'l1',
                solver = 'liblinear'))])
        clf.fit(X[train], y[train])
        results['auc'].append(roc_auc_score(y[test], clf.predict(X[test])))
        results['coef'].append(clf.named_steps['lgr'].coef_)
    return results

data_columns = np.genfromtxt('genus/text_files_for_indexing/170_columns.txt', dtype=str)
covar_encode = np.genfromtxt('genus/text_files_for_indexing/covars_ecn.txt', dtype=str)
covar_noencode = np.genfromtxt('genus/text_files_for_indexing/covars_no_ecn.txt', dtype=str)
id_columns = np.genfromtxt('genus/text_files_for_indexing/id_variables.txt', dtype=str)

# load and index data
data = pd.read_csv('brain/GENUS_FS_ATLAS_D.csv', low_memory = False)
data_use = data[data_columns]
covar_use_encode = data[covar_encode]
covar_use_noencode = data[covar_noencode]
id_use = data[id_columns]
combined = pd.concat([data_use, covar_use_encode, covar_use_noencode, id_use['IID']], axis=1)
combined = combined.dropna().drop_duplicates('IID')
brain = reidx(combined[data_columns])
cvar1 = reidx(combined[covar_encode])
cvar2 = reidx(combined[covar_noencode])
response = data[['GROUP','IID']]
response = response.set_index('IID', 1).loc[combined.IID].dropna()
response['IID'] = response.index.values
response = response.drop_duplicates('IID')
temp_encode = []
for col in cvar1.columns:
    temp_encode.append(encode(cvar1[col]))
encoded = reidx(remove_redudant(pd.concat(temp_encode, axis = 1)))
encoded = pd.concat([encoded, cvar2], axis = 1)

def y(g):
    return np.array([1. if i == 'Control' else 0. for i in g])

response = y(response['GROUP'].values)
model_data = pd.concat([brain, encoded], axis = 1)
out = logistic(X = model_data.values, y = response, splits = 7)
