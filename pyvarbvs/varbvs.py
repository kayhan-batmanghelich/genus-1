from __future__ import division
import numpy as np
import misc

def m_log_like(h):
    c = np.max(h)
    intm = np.log(np.mean(np.exp(h - c)))
    return c + intm

def bayesfactor(logw0, logw1):
    return np.exp(m_log_like(logw1) - m_log_like(logw0))

def varbvspve(X, fit, nr=1000):
    p = X.shape[0]
    ns = len(f['logw'])
    pve = np.zeros(nr)
    w = misc.normalizelogweights(fit['logw'])
    # for now a for loop...
    for i in range(nr):
        j = np.random.choice(ns, size=1, p=w)
        b = fit['mu'][:, j] + np.sqrt(fit['s'][:, j]) * np.random.standard_normal(size=p)
        b = b * (np.random.random(size=p) < fit['alpha'][:, j])
        sz = misc.var1(np.dot(X, b))
        pve[i] = sz / (sz + fit['sigma'][j])
    return pve

def varbvsnormupdate(X, sigma, sa, logodds, xy, d, alpha0, Xr0, i):
    n, p = X.shape
    for j in i:
        
