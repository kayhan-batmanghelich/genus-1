# GENUS Exploratory Analysis Repo
## This repo only reflects two primary analysis, predicting diagnosis and variational bayes

* The main results for predicting diagnosis are under the fs_cog directory. All results used an intersection of the freesurfer phenotypes and the 3 domain cognitive score(SOP, RPS, VLM) data. 9 "models" were considered for 2 classifiers - vanilla logistic regression with a sparsity inducing penalty (this is the linear classifier) and a fuzzy-trees classifier (non-linear classifier). The classifiers were chosen for the ability to 1.) reduce the dimension of the input data thus creating "interpretable results" and 2.) predict using the reduced matrix.

The models considered are:

  |Model Variable| Features
  |--------------|:--------------
  |XB            | Freesurfer phenotypes
  |XC            | Cognitive scores
  |XBC           | Freesurfer phenotypes & cognitive scores
  |XBA           | Freesurfer phenotypes & cognitive scores & covariates
  |XBCOV         | Freesurfer phenotypes & covariates
  |XCC           | Cognitive ccores & covariates
  |XBCR          | Freesurfer phenotypes with covariates projected out
  |XCCR          | Cognitive Scores with covariates projected out 


### How the classification analysis are carried out
In the <b>fs_cog/pred_diag</b> directory there is a <b>submit.py</b> file that uses two classes one for logistic regression and another for fuzzy/extra trees. These live in the <b>Mods.py</b> file which in the <b>custom utils</b> in the anaconda environment created for this project. The two classifiers inherit from a base model that simply does some minimal checks on the data, here they are: 

```python
class BaseModel(object):
    def __init__(self, X, y, out_dir, out_name):
        self.X = self.check_features(X)
        self.y = self.check_response(y)
        self.sample_check = self.check_n_samples(X, y)
        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]
        self.out = self.check_path(out_dir, out_name)

    def check_path(self, path, name):
        """
        path::string: a valid path on the file system
        name::string: name of the output file by fit_model
        checks if the path is valid, if so returns the path
        otherwise returns false.
        """
        if os.path.isdir(path):
            return os.path.join(path, name)
        else:
            error_msg = ("Specificed input path",
                         "is not a valid path")
            raise Exception(error_msg)

    def check_features(self, data):
        """
        data::pd.DataFrame: input data, the 'features'
        if input is not a pandas dataframe attemps
        to convert it to one and returns result
        """
        if not isinstance(data, pd.DataFrame):
            try:
                return pd.Dataframe(data)
            except:
                error_msg = ("X is not a pd.DataFrame and "
                             "could not convert it to a pd.DataFrame")
                raise Exception(error_msg)
        else:
            return data

    def check_response(self, response):
        """
        reponse::np.ndarray: 1d numpy array
        of the response variable
        """
        if not isinstance(response, np.ndarray):
            try:
                return np.array(response)
            except:
                error_msg = ("y is not  an np.array and could "
                             "not convert it to an anp.array")
                raise Exception(error_msg)
        else:
            return response

    def check_n_samples(self, X, y):
        """
        X::pd.DataFrame: Input data
        y::np.ndarray: response for input data
        """
        if X.shape[0] != y.shape[0]:
            error_msg = ("Features and response variable do "
                         "not have the same number of samples")
            raise Exception(error_msg)
        else:
            return True

    def save_pickle(self, data, save):
        """
        data::dict: dictionary of the object from fit_model
        save::string: path and file name of where to save output
        """
        with open(save, "w") as data_dump:
            pickle.dump(data, data_dump)
        return save

    def shuffle(self, X, y):
        """
        X::pd.DataFrame: input target values
        y::np.ndarray: response variable
        """
        X = X.values
        XY = np.hstack((X, y[:, None]))
        return XY[:, :-1], XY[:, -1]

class Logistic(BaseModel):
    def fit_model(self, X, y):
        """
        X::pd.DataFrame: Input data
        y::np.ndarray: response for input data
        """
        cv_out = StratifiedShuffleSplit(n_splits=400)
        cv_in = StratifiedKFold(n_splits=5)
        clf = Pipeline([('scaler', StandardScaler()),
                        ('lg', linear_model.LogisticRegressionCV(
                                  penalty='l1',
                                  solver='liblinear',
                                  cv=cv_in))])

        self.res = {'coef':[], 'auc':[], 'model':0}

        for idx, (train, test) in enumerate(cv_out.split(X, y)):
            clf.fit(X[train], y[train])
            prediction = clf.predict(X[test])
            self.res['coef'].append((idx, clf.named_steps['lg'].coef_[0]))
            self.res['auc'].append((idx, roc_auc_score(y[test], prediction)))

        self.res['model'] = clf
        output_saved = self.save_pickle(self.res, self.out)
        return output_saved

    def run(self):
        X, y = self.shuffle(self.X, self.y)
        return self.fit_model(X, y)

class FuzzyTrees(BaseModel):
    def fit_model(self, X, y):
        """
        X::pd.DataFrame: Input data
        y::np.ndarray: response for input data
        """
        cv_out = StratifiedShuffleSplit(n_splits=400)

        clf = Pipeline([('scaler', StandardScaler()),
                        ('fs', CustFsNoiseWinnow()),
                        ('et', ExtraTreesClassifier(n_estimators=2000))])

        self.res = {'mask':[], 'fimp':[], 'auc':[], 'model':0}

        for idx, (train, test) in enumerate(cv_out.split(X, y)):
            clf.fit(X[train], y[train])
            prediction = clf.predict(X[test])
            self.res['mask'].append((idx, clf.named_steps['fs'].mask_))
            self.res['fimp'].append((idx, clf.named_steps['et'].feature_importances_))
            self.res['auc'].append((idx, roc_auc_score(y[test], prediction)))

        self.res['model'] = clf
        output_saved = self.save_pickle(self.res, self.out)
        return output_saved

    def run(self):
        X, y = self.shuffle(self.X, self.y)
        return self.fit_model(X, y)
```
