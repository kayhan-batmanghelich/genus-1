import os
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score

data = pd.read_csv('GENUS_FS_ATLAS_D.csv', low_memory=False)
tcols = np.genfromtxt('170_columns.txt', dtype=str)
cvars = np.genfromtxt('covars_no_ecn.txt', dtype=str)
cvars = np.concatenate([cvars, ['SEX','STUDY','IID','GROUP']])
cq = ['EstimatedTotalIntraCranialVol', 'AGE_MRI', 'SEX']
sites = ['CAMH', 'CIDAR_VA', 'IMH_SIGNRP', 'UMCU_SZ2']

def res(Y):
    clinical = ('Schizophrenia', 'Clinical High Risk',
               'Familial High Risk', 'Major Depressive Disorder')
    return [1 if i in clinical else 0 for i in Y]

def ready_data(data):
    datavars = np.concatenate([tcols, cq])
    y = np.array(res(data.GROUP))
    return data[datavars], y

loso_data = {}

groups = brain_data.groupby('STUDY')

for name, group in groups:
    if name in sites:
        loso_data[name] = group

loso_dat = pd.concat(loso_data.values(), axis=0)
X, y = ready_data(loso_dat)

_, group_mask = np.unique(loso_dat.STUDY, return_inverse=True)
logo = LeaveOneGroupOut().split(X, y, group_mask)
inner_cv = StratifiedShuffleSplit(n_splits=4)
clf = Pipeline([
    ('scale', StandardScaler()),
    ('logistic', linear_model.LogisticRegressionCV(
        cv=inner_cv,
        penalty='l1',
        solver='liblinear'))
])

auc_scorer = make_scorer(roc_auc_score)
scores = cross_val_score(clf, X=X, y=y, cv=logo, scoring=auc_scorer)
