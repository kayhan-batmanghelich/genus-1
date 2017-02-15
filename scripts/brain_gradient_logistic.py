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

# make sure ids are matching and drop nans and duplicates
covar_id1 = pd.concat([id_use[['IID']],
                       covar_use_noencode],
                       axis=1).dropna().drop_duplicates('IID')
covar_id2 = pd.concat([id_use[['IID']],
                       covar_use_encode],
                       axis=1).dropna().drop_duplicates('IID')
cvars = pd.concat([covar_id1, covar_id2], axis=1).dropna()
cvars1 = cvars[covar_noencode]
cvars2 = cvars[covar_encode]
data_id = pd.concat([id_use[['IID']], data_use], axis=1)
data_id = data_id.set_index('IID', 1).loc[cvars.IID.ix[:, 0]]
data_id['IID'] = data_id.index.values
data_id = data_id.drop_duplicates(subset='IID')
id_check = np.all(data_id.IID.values == cvars.IID.ix[:, 0].values)
if not id_check:
    raise Exception("IDs not the same, cannot continue")
data_id = reidx(data_id.drop('IID', 1))
cvars1 = reidx(cvars1)
cvars2 = reidx(cvars2)
temp_encode = []
# first covariate is EstimatedTotalIntraCranialVol
for col in cvars2.columns:
    temp_encode.append(encode(cvars2[col]))
encoded_covar = pd.concat(temp_encode, axis=1)
encoded_covar = remove_redudant(encoded_covar)
covars = pd.concat([cvars1, encoded_covar], axis=1)
