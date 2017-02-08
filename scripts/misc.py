import numpy as np
import pandas as pd
from collections import Counter

def pdview(true, pred):
    t = [x for i in true for x in i]
    pred = [x for i in pred for x in i]
    df = pd.DataFrame({'y':t, 'predicted':pred})
    return df
    
def overlap(splits):
    items = (item for split in splits for item in set(split))
    return np.array([k for k, v in Counter(items).items() if v == len(splits)])

def features(names, weights):
    return map(lambda string, number: (string, round(number,5)), names, weights)

def ctype(data):
    def t(d):
        dtypes = (int, float, np.float64, np.float32, np.int64, np.int32)
        return isinstance(d, dtypes)
    try:
        iter(data)
        if np.all([t(i) for i in data]):
            return False
        else:
            return True
    except TypeError:
        if t(data):
            return False
        else:
            return True

