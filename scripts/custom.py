from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import numpy as np
import pandas as pd

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
    return np.where(~maske, avg, data)

def demean_scale(data):
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    copy = data.copy()
    copy = copy - copy.mean(0)
    return copy / copy.std(0)
