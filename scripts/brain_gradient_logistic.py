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
covar_columns = np.genfromtxt('genus/text_files_for_indexing/covariate_variables.txt', dtype=str)
id_columns = np.genfromtxt('genus/text_files_for_indexing/id_variables.txt', dtype=str)

# load and index data
data = pd.read_csv('brain/GENUS_FS_ATLAS_D.csv', low_memory = False)
data_use = data[data_columns]
covar_use = data[covar_columns]
id_use = data[id_columns]

# make sure ids are matching and drop nans and duplicates
covar_id = pd.concat([id_use[['IID']], covar_use], axis=1).dropna().drop_duplicates('IID')
data_id = pd.concat([id_use[['IID']], data_use], axis=1)
data_id = data_id.set_index('IID', 1).loc[covar_id.IID]
data_id['IID'] = data_id.index.values
data_id = data_id.drop_duplicates(subset='IID')
id_check = np.all(data_id.IID.values == covar_id.IID.values)

if not id_check:
    raise Exception("IDs not the same, cannot continue")

data_id = reidx(data_id.drop('IID', 1))
covar_id = reidx(covar_id.drop('IID', 1))
temp_encode = []
# first covariate is EstimatedTotalIntraCranialVol
for col in covar_id.columns[1:]:
    temp_encode.append(encode(covar_id[col]))
#temp_encode.append(covar_id.ix[:,0])
encoded_covar = pd.concat(temp_encode, axis=1)
