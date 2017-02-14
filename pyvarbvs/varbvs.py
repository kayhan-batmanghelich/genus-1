from __future__ import division
import numpy as np

# bayesfactor.m start
def m_log_like(h):
    c = np.max(h)
    intm = np.log(np.mean(np.exp(h - c)))
    return c + intm

def bayesfactor(logw0, logw1):
    return np.exp(m_log_like(logw1) - m_log_like(logw0))
# bayesfactor.m end

# betavar.m start
def betavar(p, mu, s):
    return p * (s + (1 - p) * mu ** 2)
# betavar.m end

# cred.m start
def cred(x, x0, w=np.ones(n)/n, c=0.95):
    n = np.shape(x)[0]
    w = w / np.sum(w)
    i, x = np.argsort(x), np.sort(x)
    # something happens to w, not sure yet
    w = w[i]
    a = np.tile(np.arange(1, n), n).reshape(n, n)
    b = a.T


# cred.m end

# diagsq.m start
def diagsq(X, a):
    return X.T.dot(a).dot(X)
# diagsq.m end

# diagsq2.m start
def diagsq2(X, A):
    return (X.dot(A) * X).sum(1)
# diagsq2.m end

# logpexp.m start
def logpexp(X):
    y = x
    i = x[x < 16]

# logpexp.m end

# int_gamma.m start
def int_gamma(logodds, alpha):
    return np.sum( (alpha - 1) * logodds + )


# normalizeweights.m start
def normalizeweights(logw):
    w = np.exp(logw - np.max(logw))
    return w / np.sum(w)
# normalizeweights.m end

# qnorm.m starts
def qnorm(x, A):
    return np.sqrt(x.dot(A).dot(x))
# qnorm.m ends
