import numpy as np
import pandas as pd

class Summary(object):
    def __init__(self, data, cols):
        self.data = data
        self.cols = cols

    def amount(self):
        locs = []
        cols = {}
        for i in self.data.columns:
            if isinstance(i, (str, np.str)):
                cols[i.lower()] = i
            elif isinstance(i, (int, np.int, float, np.float64, np.float32)):
                cols[i] = i
        for idx, val in enumerate(self.data.columns):
            if isinstance(i, (str, np.str)):
                if val.lower() in self.cols:
                    locs.append((val, idx))
        slice = self.data[[i[1] for i in locs]]
        li = [i for i in slice.values.flatten() if not pd.isnull(i)]
        rd = {i:li.count(i) for i in li}
        return rd
