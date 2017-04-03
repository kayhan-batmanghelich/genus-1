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
               verbose=True, update_sigma=True, update_sa=True, n0=0, sa0=0,
               outer_iter=np.nan):
    n, p = X.shape
    xy = np.dot(y, X)
    d = misc.diagsq(X)
    Xr = np.dot(X, (alpha*mu))
    s = sa*sigma/(sa*d + 1)
    logw = np.zeros(maxiter)
    err = np.zeros(maxiter)
    for itera in range(maxiter):
        alpha0 = alpha
        mu0 = mu
        s0 = s
        sigma0 = sigma

        logw0 = misc.int_linear(Xr, d, y, sigma, alpha, mu, s) + \
                misc.int_gamma(logodds, alpha) + \
                misc.int_klbeta(alpha, mu, s, simga*sa)

        if itera % 2:
            i = range(p)
        else:
            i = range(p, 1)
        out = varbvsnormupate(X, sigma, sa, logodds, xy, d, alpha, mu, Xr, i)
        alpha = out['alpha']
        mu = out['mu']
        Xr = out['Xr']


        logw[itera] = misc.int_linear(Xr, d, y, sigma, alpha, mu, s) + \
                     misc.int_gamma(logodds, alpha) + \
                     misc.int_klbeta(alpha, mu, s, sigma*sa)

        if update_sigma:
            sigma = (np.linalg.norm(y - Xr,ord=2)**2 +
            np.dot(d, misc.betavar(alpha, mu, s)) +
            np.dot(alpha, (s + mu**2)/sa))/(n + np.sum(alpha))

            s = sa * sigma / (sa*d + 1)

        if update_sa:
            sa = (sa0*n0 + np.dot(alpha, s + mu**2))/(n0 + sigma*np.sum(alpha))
            s = sa * sigma / (sa*d +1)

        err[itera] = np.max(np.abs(alpha - alpha0))
        if verbose:
            if np.isnan(outer_iter):
                status = np.nan
            else:
                stats = print("{}".format(outer_iter))
            # some logging functions...

        if logw[itera] < logw0:
            logw[itera] = logw0
            err[itera] = 0
            sigma = sigma0
            sa = sa0
            alpha = alpha0
            mu = mu0
            s = s0
            break
        elif err[iter] < tol:
            break

    return {'logw':logw[:itera], 'err':err[:itera],
            'sigma':sigma, 'sa':sa, 'alpha':alpha,
            'mu':mu, 's':s}


def varbvsnormindep(X, y, sigma, sa, logodds):
    s = s*sigma / (s*misc.diagsq(X) + 1)
    mu = s*np.dot(y, X) / sigma
    alpha = misc.sigmoid(logodds + (np.log(s/(sa*sigma)) + mu**2/s)/2)
    return {"alpha":alpha, "mu":mu, "s":s}

def varbvsindep(fit, X, Z, y):
    n, p = X.shape
    ns = len(fit['logw'])
    # some data shape checks...
    if np.isnan(Z).mean() == 1:
        Z = np.zeros((1, N))
        Z = np.hstack((np.ones(N)[:, None], Z))
    else:
        Z = np.hstack((np.ones(N)[:, None], Z))
    if np.any(np.isnan(y)):
        raise Exception("y cannot have any missing values")
    if fit['family'] == "gaussian":
        if Z.shape[1] == 1:
            X = X - misc.rep_row(X.mean(0), N)
            y = y - y.mean()
        else:
            #y = y - Z * (np.dot(Z.T, Z) / (Z.T * y))
            y = y - np.dot(Z, (np.dot(Z.T, Z) / np.dot(y, Z)))
            X = X - np.dot(Z, (np.dot(Z.T, Z) / np.dot(Z.T, X)))

    alpha = np.zeros((p, ns))
    mu = np.zeros((p, ns))
    s = np.zeros((p, ns))

    for i in range(ns):
        if fit["family"] == "gaussian":
            # there is a "with" statment here in the R code...
            # i think i might have to do fit["X"] instead of X and os on
            out = varbvsnormindep(X, y, sigma[i], sa[i], np.log(10)*logodds[:, i])
        elif fit["family"] == "binomial":
            out = varbvsnormindep(X, Z, y eta[, :], sa[i], np.log(10)*logoods[:, i])
        alpha[:, i] = out["alpha"]
        mu[:, i] = out["mu"]
        s[:, i] = out["s"]

    return {"alpha": alpha, "mu":mu, "s":s}
