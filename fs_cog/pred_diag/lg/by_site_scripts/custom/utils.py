import os
import numpy as np
import pandas as pd
from collections import Counter
#import h5py
import cPickle as pickle
#from pysnptools.snpreader import Bed
#import matplotlib.pyplot as plt
#from matplotlib import cm as cm
'''
def cormap(data, fs, size, cmapval, cs, title):
    """
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
    data = data.corr()
    cmap = cm.get_cmap(cs, cmapval)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(data.columns, fontsize=fs, rotation=90)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(data.index, fontsize=fs)
    cax = ax.imshow(data, interpolation="nearest", cmap=cmap, clim=(-1,1))
    fig.colorbar(cax, fraction=0.033, pad=0.035)
    fig.set_size_inches(*size)
    plt.title(title, size=20)
    return None
'''
domain_scores = ('SOPdomainAvgZ', 'ATVIdomainAvgZ', 'VWMdomainAvgZ',
        'NVWMdomainAvgZ', 'VLMdomainAvgZ', 'NVLMdomainAvgZ',
        'RPSdomainAvgZ', 'VISPAdomainAvgZ')

data_values = ('cog', 'gen', 'fam')

abbrev_1 = """camh cid1 cid2 deco gapl imhs lanr mci1
              mci2 mci3 mfs1 mts1 nef1 nef2 nef3 nef4
              nef5 phrs tcin nuig usz1 usz2""".split()

abbrev_2 = """camh cati cida cogs deco gapl imhs imh4
              lanr mcic mfs1 mts1 nefs paco page phrs
              tcdn umcu""".split()

def id_intersect(data):
    ids = set(data[0]).intersection(*data)
    return np.array(list(ids))

def gene_string_check(path):
    t = (".bed", ".fam", ".bim")
    m = path[-4:]
    if (t[0] == m or t[1] == m or t[2] == m):
        path = path[:-4]
    return path

def ids_and_fam(path):
    path = gene_string_check(path)
    bed_IDs = Bed(path + ".bed").iid.astype(str)
    print("Bed::\nCol 1 FID\nCol 2 IID\n")
    fam_data = pd.read_csv(path + ".fam", sep=' ', header=None)
    fam_data.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'DIAG']
    fam_data['IID'] = fam_data['IID'].astype(str)
    print("DIAG::\n0: Unknown\n1: Control\n2: Case\n")
    print("SEX::\n1: Male\n2: Female")
    assert((bed_IDs[:, 1] == fam_data.IID.values).mean()== 1.)
    return fam_data, bed_IDs

def pysnp_genpreproc(data):
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    nan_locs = np.isnan(data)
    data[nan_locs] = 5
    mask = data != 5
    avg = np.true_divide((data*mask).sum(0), mask.sum(0))
    return np.where(~mask, avg, data)

def h5_save(path, data_obj, dt):
    with h5py.File(path, "w") as file_store:
        for key, val in data_obj.items():
            data_set = file_store.create_dataset(
                key, val.shape, dtype=dt)
            data_set[...] = val
    return path, data_obj.keys()

def h5_read(path, key):
    with h5py.File(path, "r") as file_store:
        data = file_store[key][...]
    return data

def read_pickle(path):
    with open(path, "r") as file_store:
        data = pickle.load(file_store)
    return data

def save_pickle(path, data):
    with open(path, "w") as file_store:
        pickle.dump(data, file_store)
    return path

def make_non_singular(X, tol = 1e-05):
    Q, R = np.linalg.qr(X)
    independent = np.where(np.abs(R.diagonal()) > tol)[0]
    return X[:, independent]

def proj(X , C):
    P = np.eye(C.shape[0]) - C.dot(np.linalg.inv(C.T.dot(C))).dot(C.T)
    return P.dot(X)

def numeric_func_inter(x):
    try:
        x.astype(float)
        return True
    except:
        return False

def get_numeric_col_only(X):
    if isinstance(X, pd.DataFrame): 
        X = X.values
    return map(numeric_func_inter,
              [X[:, i] for i in range(X.shape[1])]) 

def hdfload(path, dtype, key):
    print("Returning train, test...")
    train = pd.read_hdf(path, '{}_train_{}'.format(dtype, key))
    test = pd.read_hdf(path, '{}_test_{}'.format(dtype, key))
    return train, test

class Summary(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def fit(self):
        to_count = Counter(np.array([self.data[i].values for i in 
                             self.columns]).flatten().tolist())
        return {i: to_count[i] for i in to_count if not pd.isnull(i)}

def encoder(X):
    n = len(X)
    uniq, inv = np.unique(X, return_inverse=True)
    result = np.zeros((n, len(uniq)), dtype=float)
    result[np.arange(n), inv] = 1
    return result

def remove_redudant(encoded):
    to_drop = []
    for col in range(1, encoded.shape[1]):
        if (encoded[:, col - 1] + encoded[:, col]).mean() == 1:
            to_drop.append(False)
        else:
            to_drop.append(True)
    return encoded[:, to_drop]

def get_nonzerocoef_cols(cols, dobj, idx_name, idx_num):
    return cols[np.nonzero(dobj[idx_name][idx_num])[1]]

def genomic_preproc(data):
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    mask = data != 3
    avg = np.true_divide((data * mask).sum(0), mask.sum(0))
    return np.where(~mask, avg, data)

def demean_scale(data):
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    copy = data.copy()
    copy = copy - copy.mean(0)
    return copy / copy.std(0)    
    
class Match(object):
    def __init__(self, id_var, data):
        self.data = data
        self.id_var = id_var

    def check_input(self, input_list):
        first_pass = lambda x: True if isinstance(x, pd.DataFrame) else False
        return [item for item in input_list if first_pass(item)]

    def dropnans(self):
        data_nonans = []
        for data_type in self.data:
            data_nonans.append(data_type.dropna(0))
        return data_nonans

    def inter(self, loi):
        return list(set(loi[0]).intersection(*loi))

    def get_matching_ids(self, data):
        data = self.check_input(data)
        return self.inter([d[self.id_var].values for d in data ])

    def index(self, id_var, data, ids):
        data = data.set_index(id_var).loc[ids].reset_index()
        return data.drop_duplicates(subset=id_var).dropna()
        
    def fit(self):
        new_data = self.dropnans()
        ids = self.get_matching_ids(new_data)
        out = []
        for data in new_data:
            out.append(self.index(self.id_var, data, ids))
        return out

def process_data(input_path, output_path, name):
    snpreader = Bed(os.path.join(input_path, name))
    data = snpreader.read()
    values = data.val
    preproc_vals = pysnp_genpreproc(values)
    assert(np.any(np.isnan(preproc_vals)) == False)
    saved = os.path.join(output_path, name + ".h5py")
    path, keys = h5_save(path=saved, data_obj={name:preproc_vals}, dt='f')
    return {'n_subjects':data.iid_count, 'subject_ids':data.iid,
            'n_snps':data.sid_count, 'snp_ids':data.sid,
            'data_preprocessed_location': {'path':path, 'key':keys}}
