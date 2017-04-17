import os
import numpy as np
import nibabel as nb
from surfer import Brain

ldir = '/cm/shared/openmind/freesurfer/5.3.0/subjects/fsaverage/label'
lfile = 'lh.aparc.a2009s.annot'
aparc_file = os.path.join(ldir, lfile)
labels, ctab, names = nb.freesurfer.read_annot(aparc_file)

#brain = Brain('fsaverage', 'lh', 'inflated', background="white")


brain = Brain("fsaverage", "split", "inflated",
              views=['lat', 'med'], background="white")


import matplotlib.pyplot as plt
import pickle
from collections import Counter

def read_pickle(name):
    with open(name, "rb") as data:
        data = pickle.load(data)
    return data

results = '/storage/gablab001/data/genus/pkls/rf_results.pkl'
res = read_pickle(results)
#txtb = '/storage/gablab001/data/genus/current/structured/genus/text_files_for_indexing'
#bd = '/storage/gablab001/data/genus/current/structured/brain/'
#thickness = np.genfromtxt(os.path.join(txtb, '170_columns.txt'), dtype=str)

def make_dataframe(x):
    dflist = []
    cols=['iteration', 'feat_index', 'hemi','feat', 'auc']
    for key, val in x.items():
        split_key = key.split('_')
        if 'lh' or 'rh' not in split_key:
            if 'rh' in split_key[2]:
                dflist.append([int(split_key[0]), int(split_key[1]),
                           'rh', '_'.join(split_key[2:]), val])
            elif 'lh' in split_key[2]:
                dflist.append([int(split_key[0]), int(split_key[1]),
                           'lh', '_'.join(split_key[2:]), val])
        else:
            dflist.append([int(split_key[0]), int(split_key[1]),
                        split_key[2], '_'.join(split_key[2:]), val])
    return pd.DataFrame(dflist, columns=cols)

results_df = make_dataframe(res)
groups = results_df.groupby('iteration')
smallest = {'rh':[[],[]], 'lh': [[],[]]}
for i in range(100):
    itera = groups.get_group(i)
    hemis = itera.groupby('hemi')
    for hemi in ['lh','rh']:
        data = hemis.get_group(hemi)
        gir = data[data['auc'] == data['auc'].min()]
        smallest[hemi][0].append(gir['feat'].values[0].replace('lh_','') \
                                 .replace('rh_','') \
                                 .replace('_thickness_D', ''))
        smallest[hemi][1].append(gir['auc'].values[0])

rh_count = dict(Counter(smallest['rh'][0]))
right_df = pd.DataFrame({'col':[i.replace('.', '-') for i in rh_count.keys()],
                         'val': rh_count.values()})
right_df = right_df.set_index('col').loc[names].reset_index().fillna(0)
vtx_rh = right_df.val.values[labels]

lh_count = dict(Counter(smallest['lh'][0]))
left_df = pd.DataFrame({'col':[i.replace('.', '-') for i in lh_count.keys()],
                         'val': lh_count.values()})
left_df = left_df.set_index('col').loc[names].reset_index().fillna(0)
vtx_lh = left_df.val.values[labels]

brain.add_data(vtx_rh,
               right_df.val.min(),
               right_df.val.max(),
               colormap="Blues",
               alpha=.8,
               hemi='rh')

brain.add_data(vtx_lh,
               left_df.val.min(),
               left_df.val.max(),
               colormap="Blues",
               alpha=.8,
               hemi='lh')

brain.add_annotation(aparc_file, hemi='lh')
brain.add_annotation(aparc_file, hemi='rh', remove_existing=False)
