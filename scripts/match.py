import numpy as np
import pandas as pd
import magic

class Match(object):
    def __init__(self, id_var, brain=[], cognitive=[], genomic=[],  sep = None):
        self.brain = brain
        self.cognitive = cognitive
        self.genomic = genomic
        self.id_var = id_var
        self.sep = sep

    def inter(self, loi):
        return list(set(loi[0]).intersection(*loi))
    
    def has_item(self, x):
        try:
            x[0] 
        except IndexError:
            return False
        return True
    
    def check_type(self, x):
        if isinstance(x, pd.DataFrame):
            self.has_item(x)
        else:
            return False
    
    def get_matching_ids(self, id_var, b, c, g):
        
        id_var_inter = self.inter([val.columns.values for val in (b, c, g) \
                                  if self.check_type(val.columns.values)])
        if not id_var_inter[0] == id_var:
            raise Exception("Some input data is missing ID variable")
        else:
            reduce_n = self.inter([val[id_var_inter[0]] for val in (b, c, g) \
                                  if self.check_type(val[id_var_inter[0]].values)])
            return reduce_n

    def load(self, data):
        try:
            dl = magic.from_file(data)
            print dl
            if 'Hierarchical Data Format' in dl:
                return pd.read_hdf(data)
            elif 'ASCII' in dl:
                if not self.sep:
                    return pd.read_csv(data, sep = self.sep)
                return pd.read_csv(data)
            elif 'Matlab' or 'mat-file' in dl:
                return scipy.io.loadmat(data)
        except TypeError:
            pass

    def index(self, id_var, data, ids):
        return data.set_index(id_var, 1).loc[ids]

    def out(self):
        data = {
        'b': self.load(self.brain),
        'c': self.load(self.cognitive),
        'g': self.load(self.genomic)
        }

        ids = self.get_matching_ids(self.id_var, **data)

        return (self.index(self.id_var, data['b'], ids),
                self.index(self.id_var, data['c'], ids),
                self.index(self.id_var, data['g'], ids))
