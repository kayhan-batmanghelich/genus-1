import numpy as np
import pandas as pd
from collections import Counter

class Summary(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def fit(self):
        to_count = Counter(np.array([self.data[i].values for i in self.columns]).flatten().tolist())
        return {i: to_count[i] for i in to_count if not pd.isnull(i)}

def encode(data):
    encoded_label = data.name
    le = LabelEncoder()
    ohe = OneHotEncoder()
    temp_encode = le.fit_transform(data)
    final_encode = ohe.fit_transform(temp_encode.reshape(-1,1)).toarray()
    return pd.DataFrame(final_encode, columns = le.classes_)

def remove_redudant(encoded):
    to_drop = []
    for col in range(1, encoded.shape[1]):
        if (encoded.ix[:, col - 1] + encoded.ix[:, col]).mean() == 1:
            to_drop.append(encoded.columns[col])
    return encoded.drop(to_drop,axis=1)

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

    def inter(self, loi):
        return list(set(loi[0]).intersection(*loi))

    def get_matching_ids(self, id_var, data):
        data = self.check_input(data)
        id_var = self.inter([i.columns.values for i in data])
        return self.inter([d[id_var[0]].values for d in data ])

    def index(self, id_var, data, ids):
        data = data.set_index(id_var, 1).loc[ids]
        data[id_var] = data.index.values
        return data.drop_duplicates(id_var).drop(id_var, 1)
        
    def fit(self):
        ids = self.get_matching_ids(self.id_var, self.data)
        out = []
        for data in self.data:
            out.append(self.index(self.id_var, data, ids))
        return out
