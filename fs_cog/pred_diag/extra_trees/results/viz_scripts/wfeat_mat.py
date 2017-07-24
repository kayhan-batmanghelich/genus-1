import os
import numpy as np
import pandas as pd
from custom import utils
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import cm as cm

def featmap(data, title, sn, fs=20, size=(15,15), cmapval=100, cs="jet"):
    """
    data: pandas.DataFrame
    fs: font size::int
    size: size of graph::tuple of length 2
    cmapval: n discrete color breaks in cs::int
    cs: color scheme::string
    title: graph title::string
    """
    if not isinstance(data, pd.DataFrame):
        raise Exception("data must be a pandas.DataFrame")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1) 
    cmap = cm.get_cmap(cs, cmapval)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(data.columns, fontsize=fs, rotation=90)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(data.index, fontsize=fs)
    cax = ax.imshow(data, interpolation="nearest", cmap=cmap, clim=(0,.5))
    fig.colorbar(cax, fraction=0.033, pad=0.035)
    fig.set_size_inches(*size)
    plt.title(title, size=20)
    plt.savefig(sn, bbox_inches="tight", dpi=300)
    plt.close()
    return None
"""
task:
create a weighted feature importance matrix,
where the weights are the AUC score of that data point
in the distribution

weighted_importance = sum_i(auc_i*feat_import_vect)/sum_i(auc_i)

afterwards:
do argsort to put most important on the left
"""
# for later u# header path
hp = "/storage/gablab001/data/genus/fs_cog/pred_diag/column_headers"
dp = "/storage/gablab001/data/genus/fs_cog/pred_diag/extra_trees/results/"
dpv = [i for i in os.listdir(dp) if ".pkl" in i]
dnh = {item:item.split('_')[1].replace(".pkl", "") for item in dpv}

def equate_n_sort(abso, subset):
    if len(abso) < len(subset):
        raise Exception("first arg must be the full set")
    diff = np.setdiff1d(abso, subset)
    res = np.concatenate([subset, diff])
    res.sort()
    if not np.all(res == sorted(abso)):
        raise Exception("could not equate inputs")
    else:
        df = pd.DataFrame(columns=res)
        df.loc[0, :] = np.zeros(len(res))
        return df

# header for big matrix
all_headers = np.genfromtxt(os.path.join(hp, "XBA"), dtype=str)
data_save = []

#def extract(res, k2):
#    k1 = k2[-1]
#    tmp = [(int(key.split(k1)[1]), val) for 
#           key, val in res.items() if k2 in key]
#    return sorted(tmp, key=lambda tup: tup[0])

#def mask2imp(mask, imp):
#    max_mask = len(mask)
#    out_imp = []
#    imp = imp.tolist()
#    for idx in range(max_mask):
#        if mask[idx]:
#            out_imp.append(imp[idx])
#        else:
#            imp.insert(idx, 0)
#            out_imp.append(0)
#    return np.array(out_imp)

def mask2imp(mask, imp):
    mask_copy = mask.copy().astype(float)
    mask_copy[mask_copy > 0] = imp
    return mask_copy

for key_i, val_i in dnh.items():
    result = utils.read_pickle(os.path.join(dp, key_i))
    header = np.genfromtxt(os.path.join(hp, val_i), dtype=str)
    #fimp = extract(result, "fimp")
    #aucs = [i[1] for i in extract(result, "auc")]
    #masks = extract(result, "mask")
    fimp = np.array([f[1] for f in result['fimp']])
    aucs = [a[1] for a in result['auc']]
    masks = np.array([m[1] for m in result['mask']])
    fimp_mat = np.array([mask2imp(masks[i], fimp[i]) for i in range(len(fimp))])
    #fimp_mat = np.array([fimp[i] for i in range(len(fimp))])
    wf_imp =(np.abs(fimp_mat)*np.array(aucs)[:, None]).sum(0)/np.sum(aucs)
    df = pd.DataFrame(columns=header)
    df.loc[0, :] = wf_imp
    fn = equate_n_sort(all_headers, header)
    data_save.append(pd.concat([fn, df]).fillna(0).sum(0))

def most_imp(data, nth_p):
    perc = np.percentile(data.sum(0).values, nth_p)
    return data.iloc[:, (data.sum(0) > perc).values]

matrix = pd.concat(data_save, axis=1)
matrix.columns=dnh.values()
matrix = matrix.T

row_sort = ['XBA','XBC','XBCCR','XBCOV','XBCR','XB','XCC','XCCR','XC']
matrix['icol'] = matrix.index.values
matrix = matrix.set_index('icol').loc[row_sort]
sort_idx = np.argsort(matrix.sum(0).values)[::-1]
matrix = matrix.iloc[:, sort_idx]

for percentile in [95, 90, 85, 80]:
    mat = most_imp(matrix, percentile)
    title = "Importance at {} percentile".format(percentile)
    save_name = "feat_mat_et_{}.png".format(percentile)
    featmap(mat, title, save_name)

