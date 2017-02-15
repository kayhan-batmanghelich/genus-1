import numpy as np
import pandas as pd

class Summary(object):

    sex_cols = ('sex', 'gender', 'male', 'female')

    def __init__(self, data):
        self.data = data

    def loc(self):
        locs = []
        cols = {i.lower():i for i in self.data.columns}
        for idx, val in enumerate(self.data.columns):
            if val.lower() in self.sex_cols:
                locs.append((val, idx))
        return locs

    def intermediate(self, locs):
        locs = [i[1] for i in locs]
        slice = self.data[locs]
        li = [i for i in slice.values.flatten() if not pd.isnull(i)]
        num_sex = len(set(li))
        nan = 0
        for i in slice.values.flatten():
            if pd.isnull(i):
                nan += 1
        rd = {i:li.count(i) for i in li}
        rd['nan'] = nan
        rd['nsex'] = num_sex
        return rd
