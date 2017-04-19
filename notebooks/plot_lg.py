import os
import numpy as np
import nibabel as nb
from surfer import Brain

ldir = '/cm/shared/openmind/freesurfer/5.3.0/subjects/fsaverage/label'
lfile = 'lh.aparc.a2009s.annot'
aparc_file = os.path.join(ldir, lfile)
labels, ctab, names = nb.freesurfer.read_annot(aparc_file)

brain = Brain("fsaverage", "split", "inflated",
              views=['lat', 'med'], background="white")

import matplotlib.pyplot as plt
import pickle
from collections import Counter
import pandas as pd

def read_pickle(name):
    with open(name, "rb") as data:
        data = pickle.load(data)
    return data

left = read_pickle('/storage/gablab001/data/genus/pkls/lh_100_iters_nonzero.pkl')
right = read_pickle('/storage/gablab001/data/genus/pkls/rh_100_iters_nonzero.pkl')
left_df = pd.DataFrame({'col':left.keys(), 'val': left.values()})
right_df = pd.DataFrame({'col':right.keys(), 'val': right.values()})
left_df = left_df.set_index('col').loc[names].reset_index()
right_df = right_df.set_index('col').loc[names].reset_index()

vtx_data = left_df.val.values[labels]
vtx_rh = right_df.val.values[labels]

brain.add_data(vtx_data,
               left_df.val.min(),
               left_df.val.max(),
               colormap="Reds",
               alpha=.8,
               hemi='lh')

brain.add_annotation(aparc_file, hemi='lh')

brain.add_data(vtx_rh,
               right_df.val.min(),
               right_df.val.max(),
               colormap="Reds",
               alpha=.8,
               hemi='rh')

brain.add_annotation(aparc_file, hemi='rh', remove_existing=False)
