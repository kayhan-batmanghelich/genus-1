from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

def encode(data):
    encoded_label = data.name
    le = LabelEncoder()
    ohe = OneHotEncoder()
    temp_encode = le.fit_transform(data)
    final_encode = ohe.fit_transform(temp_encode.reshape(-1,1)).toarray()
    return pd.DataFrame(final_encode, columns = le.classes_)

def res(y):
    p = 'Schizophrenia'
    return np.array([1. if val == p else 0 for val in y])

def remove_redudant(encoded):
    to_drop = []
    for col in range(1, encoded.shape[1]):
        if (encoded.ix[:, col - 1] + encoded.ix[:, col]).mean() == 1:
            to_drop.append(encoded.columns[col])
    return encoded.drop(to_drop,axis=1)

def reidx(df):
    df1 = df.copy()
    df1.index = [x for x in range(df.shape[0])]
    return df1
