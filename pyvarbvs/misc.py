from __future__ import division
import numpy as np

eps = np.finfo(float).eps

def tf2yn(x):
    if x:
        return True
    return False

def var1(x):
    n = len(x)
    return np.var(x) * (n - 1) / n

def var_cols(X):
    return  np.array(map(var1, X.T))

def qnorm(x, a):
    if len(a.shape) > 1:
        return np.sqrt(x.dot(a).dot(x))
    return np.sqrt((x*a).dot(x))

def rep_row(x, n):
    return np.tile(x, n).reshape(n, len(x)).T

def diagsq(X, a=None):
    if not a:
        return np.diag(X.T.dot(X))
    return np.diag(X.T.dot(np.diag(a)).dot(X))

def diagsqt(X, a=None):
    if not a:
        return np.diag(X.T.dot(X))
    return np.diag(X.dot(np.diag(a)).dot(X.T))

def diagsq2(X, A):
    return (X.dot(A) * X).sum(1)

def sigmoid(x):
    return 1. / (1 + np.exp(-x))

def logit(x):
    return np.log((x + eps) / ((1 - x) + eps))

def logpexp(x):
    y = x.copy().astype(float)
    i = np.where(x < 16)[0]
    y[i] = np.log(1 + np.exp(x[i]))
    return y

def logsigmoid(x):
    return -logpexp(-x)

def slope(x):
    return (sigmoid(x) - .5) / (x + eps)

def int_gamma(logodds, alpha):
    left = (alpha - 1) * logodds
    right = logsigmoid(logodds)
    return np.sum(left + right)

def int_klbeta(alpha, mu, s, sa):
    return (np.sum(alpha) + np.dot(alpha, np.log(s / sa)) - \
            np.dot(alpha, s + mu**2) / sa) / 2 - \
            np.dot(alpha, np.log(alpha + eps)) - \
            np.dot(1 - alpha, np.log(1 - alpha + eps))

def betavar(p, mu, s):
    return p * (s + (1 - p) * mu**2)

def normalizelogweights(logw):
    c = np.max(logw)
    w = np.exp(logw - c)
    return w / w.sum()
