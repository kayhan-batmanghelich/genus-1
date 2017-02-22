import numpy as np
import pandas as pd

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
