import numpy as np
import pandas as pd
import magic

class Match(object):
    def __init__(self, brain, genomic, cognitive, id_var, sep):
        self.brain = brain
        self.genmoic = genomic,
        self.cognitive = cognitive
        self.id_var = id_var,
        self.sep = sep

    data_store_types = ('Hierarchical Data Format',
                        'ASCII', 'Matlab', 'mat-file')
    def inter(self, loi):
        return list(set(loi[0]).intersection(*loi))

    def get_matching_ids(self, id_var, b, c, g):
        id_var_inter = inter([val.columns.values for val in (b, c, g)])
        if not id_var_inter[0] == id_var:
            raise Exception("Some input data is missing ID variable")
        else:
            reduce_n = inter([val[id_var_inter[0]] for val in (b, c, g)])
            return reduce_n

    def load(self, data, sep = ','):
        dl = magic.from_file(data)
        if 'Hierarchical Data Format' in dl:
            return pd.read_hdf(data)
        elif 'ASCII' in dl:
            return pd.read_csv(data, sep = sep)
        elif 'Matlab' or 'mat-file' in dl:
            return scipy.io.loadmat(data)

    def index(self, id_var, data, ids):
        return data.set_index(id_var).loc(ids)

    def out(self, id_var, b, c, g, ids):
        
        data = {'b': self.brain, 'c': self.cognitive, 'g': self.genomic}
        ids = get_matching_ids(self.id_var, **data)

        return (index(self.id_var, self.brain, self.id_var),
                index(self.id_var, self.cognitive, self.id_var),
                index(self.id_var, self.genomic, self.id_var))
