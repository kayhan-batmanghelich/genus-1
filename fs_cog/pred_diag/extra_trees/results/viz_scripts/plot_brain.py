import os
import sys
import numpy as np
from custom import utils
import nibabel as nb
from surfer import Brain

result_name = sys.argv[1]
header_name = sys.argv[2]
save_name = sys.argv[3]

brain = Brain("fsaverage", "split", "inflated",
               views=['lat', 'med'], background="white")


hp = "/storage/gablab001/data/genus/fs_cog/pred_diag/column_headers"
rp = "/storage/gablab001/data/genus/fs_cog/pred_diag/extra_trees/results/"
with open(os.path.join(hp, header_name), "r") as tmp:
    headers = np.array(tmp.read().strip("\n").split("\n"))

res = utils.read_pickle(os.path.join(rp,"{}".format(result_name)))

print("loading result: {}".format(result_name))
print("loading header: {}".format(header_name))

coefs = [m[1] for m in res["mask"]]
features_selected_lh = []
features_selected_rh = []
features_all = []

from collections import Counter
import pandas as pd

for coef in coefs:
    features_all.extend(headers[coef].tolist())
    for col in headers[coef].tolist():
        if 'lh' in col:
            features_selected_lh.append(col)
        elif 'rh' in col:
            features_selected_rh.append(col)

left_count = dict(Counter(features_selected_lh))
right_count = dict(Counter(features_selected_rh))
left_df = pd.DataFrame({
    'col': [key.replace('lh_','').replace('.','-').replace('_thickness_D','') 
            for key, val in left_count.items()],
    'val': [val for key, val in left_count.items()]
})
right_df = pd.DataFrame({
    'col': [key.replace('rh_','').replace('.','-').replace('_thickness_D','') 
            for key, val in right_count.items()],
    'val': [val for key, val in right_count.items()]
})
label_dir = '/cm/shared/openmind/freesurfer/5.3.0/subjects/fsaverage/label'
left_label_file = 'lh.aparc.a2009s.annot'
right_label_file = 'rh.aparc.a2009s.annot'
lh_aparc_file = os.path.join(label_dir, left_label_file)
rh_aparc_file = os.path.join(label_dir, right_label_file)
lh_labels, lh_ctab, lh_names = nb.freesurfer.read_annot(lh_aparc_file)
rh_labels, rh_ctab, rh_names = nb.freesurfer.read_annot(rh_aparc_file)
left_df = left_df.set_index('col').loc[lh_names].reset_index().fillna(0)
right_df = right_df.set_index('col').loc[rh_names].reset_index().fillna(0)
vtx_lh = left_df.val.values[lh_labels]
vtx_lh[lh_labels == -1] = 0
vtx_rh = right_df.val.values[rh_labels]
vtx_rh[rh_labels == -1] = 0

brain.add_data(vtx_lh,
               0,
               400,
               colormap="Reds",
               alpha=.8,
               hemi='lh')

brain.add_annotation(lh_aparc_file, hemi='lh')

brain.add_data(vtx_rh,
               0,
               400,
               colormap="Reds",
               alpha=.8,
               hemi='rh')

brain.add_annotation(rh_aparc_file, hemi='rh', remove_existing=False)
save_name = "../images/et_{}_brain.png".format(save_name)
brain.save_image(save_name)

