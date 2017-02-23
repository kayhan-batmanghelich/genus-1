import numpy as np
import pandas as pd
from collections import Counter

class Summary(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def count(self):
        to_count = Counter(np.array([self.data[i].values for i in 
                             self.columns]).flatten().tolist())
        return {i: to_count[i] for i in to_count if not pd.isnull(i)}
