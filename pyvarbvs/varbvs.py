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
        s = sa * sigma / (sa * d[j] + 1)
        r = alpha[j] + mu[j]
        mu[j] = s / sigma * (xy[j] + d[j]*r - np.sum(x[:,j]*Xr))
        alpha[j] = misc.sigmoid(logodds[j] + (np.log(s/(sa*sigma)) + mu[j]**2/s)/2)
        Xr = Xr + (alpha[j] * mj[j] - r) * X[:, j]
    return {'alpha':alpha, 'mu':mu, 'Xr':Xr}

def varbvsnorm(X, y, sigma, sa, logodds, alpha, mu, tol=1e-4, maxiter=1e4,
               verbose=True, update_sa=True, n0=0, sa0=0):
    n, p = X.shape
    xy = np.dot(y, X)
    d = misc.diagsq(X)
    Xr = np.dot(X, (alpha*mu))
    s = sa*sigma/(sa*d + 1)
    logw = np.zeros(maxiter)
    err = np.zeros(maxiter)
    for iter in range(maxiter):
        alpha0 = alpha
        mu0 = mu
        s0 = s
        sigma0 = sigma

        logw0 = misc.int_linear(Xr, d, y, sigma, alpha, mu, s) + \
                misc.int_gamma(logodds, alpha) + \
                misc.int_klbeta(alpha, mu, s, simga*sa)

        if iter % 2:
            i = range(p)
        else:
            i = range(p, 1)
        out = varbvsnormupate(X, sigma, sa, logodds, xy, d, alpha, mu, Xr, i)
        alpha = out['alpha']
        mu = out['mu']
        Xr = out['Xr']


        logw[iter] = misc.int_linear()
