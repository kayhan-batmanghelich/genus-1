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
```
