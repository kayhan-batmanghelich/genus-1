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
