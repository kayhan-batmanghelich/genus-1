import os
import numpy as np
import pandas as pd
from custom import utils


# VALIDATION DATA
valbase = '/storage/gablab001/data/genus/current/structured/validate'

# 183_cor_subcor_D_thickness.csv contains the 170 features for the validation set
valdatacob = '183_cor_subcor_D_thickness.csv'
valdatacob = pd.read_csv(os.path.join(valbase, valdatacob))

# the file containing the covariates
valcov = pd.read_csv(os.path.join(valbase, '183_covariates_to_use.csv'))

# this is y in the classification analysis
valres = [1 if i == 0 else 0 for i in valcov['diag'].values]

valcov = valcov[['ICV','age','gender']]

# base directory where some text files live, they incude header names
# example: example_header.txt will have:
# pallidum
# cerebellum
txtb = '/storage/gablab001/data/genus/current/structured/genus/text_files_for_indexing'

# directory to where the genus brain that is
bd = '/storage/gablab001/data/genus/current/structured/brain/'

# the header text file for the 170 columns, this will be used to subset the 
# entire genus data
thickness = np.genfromtxt(os.path.join(txtb, '170_columns.txt'), dtype=str)

# making equivalent volume headers from the thickness headers
volume = ' '.join(thickness).replace('thickness_D','volume_D').split(' ')

# this variable is for convenience so that i dont hav to change
# where the headers are in multiple places, just here
colheads = thickness

# loading the covariate headers that will be one hot encoded
cvar_encode = np.genfromtxt(os.path.join(txtb, 'covars_ecn.txt'), dtype=str)

# the covariate headers that wont be one hot encoded
cvar = np.genfromtxt(os.path.join(txtb, 'covars_no_ecn.txt'), dtype=str)

# GENUS brain data
brain = pd.read_csv(os.path.join(bd, 'GENUS_FS_ATLAS_D.csv'), low_memory=False)

# GENUS response variable
response = brain[['IID','GROUP']]

# here i combined all the needed data so that I can drop rows all together and 
# make sure all parts of the data, brain regions, covariates, ID, response are
# sorted by the same rows
combined = pd.concat([
    brain[colheads],
    brain[cvar],
    brain[cvar_encode],
    brain[['IID','GROUP']]
], axis=1).dropna().drop_duplicates('IID')

# the genus y in the classification analysis
response = combined['GROUP'].values

# subsetting genus data to only include the 170 features
X_data = combined[colheads].reset_index(drop=True)

# getting the covariates that wont be one hot encoded
cvar_ne = combined[cvar].reset_index(drop=True)

# and the covariates that will be one hot encoded
cvar_e = combined[cvar_encode].reset_index(drop=True)

# performing one hot encoding
cvar_e = pd.concat([
    pd.DataFrame(utils.encoder(cvar_e[col])) for col in cvar_e.columns
], axis=1, ignore_index=True)

# recombining covariates
cvars = pd.concat([cvar_ne, cvar_e], axis=1)

# response variable to be fed into classifier
y = np.array([1 if i == 'Schizophrenia' else 0 for i in response])

# remove effects of covariates
# NOTE: this will result in extremely poor results for the validation set
#       because the same covariates cannot be applied to the validation set
#       it's only when the same covariates are removed from each set that we see
#       results similar to what I show
non_sing_cvars = utils.make_non_singular(cvars.values)
XG = utils.proj(X_data.values, non_sing_cvars)
XCO = utils.proj(valdatacob.values, valcov.values)

# dictionary that will be passed to the classifier node
data_dict = {'GENUS_X':XG, 'GENUS_y':y,
             'GENUS_XCOLS':X_data.columns.values,
             'GENUS_COVCOLS':cvars.columns.values,
             'COBREFMRI_X': XCO, 'COBREFMRI_y': valres,
             'COBREFMRI_XCOLS': valdatacob.columns.values,
             'COBREFMRI_COVCOLS': valcov.columns.values}

from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit

cv_in = StratifiedKFold(n_splits=5, shuffle=True)
cv_out = StratifiedShuffleSplit(n_splits=400)

clf = Pipeline([
    ('scale', StandardScaler()),
    ('lg', linear_model.LogisticRegressionCV(
        penalty='l1',
        solver='liblinear',
        cv=cv_in,
        Cs=200,
        n_jobs=-1
    ))
])

results = {}

for idx, (train, test) in enumerate(cv_out.split(data_dict['GENUS_X'], data_dict['GENUS_y'])):
    
    clf.fit(data_dict['GENUS_X'][train], data_dict['GENUS_y'][train])
    
    results['coef_{}'.format(idx)] = clf.named_steps['lg'].coef_
    
    results['test_auc_split_{}'.format(idx)] = roc_auc_score(
                                             data_dict['GENUS_y'][test], 
                                             clf.predict(data_dict['GENUS_X'][test])
                                         )
    results['train_auc_split_{}'.format(idx)] = roc_auc_score(
                                              data_dict['GENUS_y'][train],
                                              clf.predict(data_dict['GENUS_X'][train])
                                         )
    
    
import pickle

with open("results_logisticregression.pkl", "wb") as f:
    pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
