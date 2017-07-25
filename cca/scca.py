import numpy as np

def partial_cov(X, Y):
    n = X.shape[0]
    cX = X - np.sum(X, axis=0, keepdims=True) / n
    cY = Y - np.sum(Y, axis=0, keepdims=True) / n
    cov_XY = 1. / (n - 1) * np.dot(cX.T, cY)
    return cov_XY 

def mean_center_scale(x):
    x = x - x.mean(0)
    return x / x.std(0)

class SCCA:
    """Implementation of sparce canonical correlation
    Reference: Parkhomenko et al.: Sparse Canonical Correlation Analysis
    This is a work in progress, not expected to work yet.
    """
    def __init__(self, X, Y, c1=0.01, c2=0.01, lu=0.00000001, lv=0.00000002):
        self.X = mean_center_scale(X)
        self.Y = mean_center_scale(Y)
        self.p = np.shape(X)[1]
        self.q = np.shape(Y)[1]
        self.c1 = c1
        self.c2 = c2
        self.lu = lu
        self.lv = lv

    def _K(self):
        xy = partial_cov(self.X, self.Y)
        xxi = np.linalg.inv(np.diag(np.diag(np.cov(X, rowvar=False))))
        yyi = np.linalg.inv(np.diag(np.diag(np.cov(Y, rowvar=False))))
        self.xinv = xxi
        self.yinv = yyi
        return np.dot(xxi, xy).dot(yyi)

    def _norm(self, w):
        lw=np.sqrt(np.dot(w, w))
        return w/lw

    def _soft_thresh(self, w, t):
        if t == 'u':
            w = (np.abs(w) - .5*self.lu) + np.sign(w)
        elif t == 'v':
            w = (np.abs(w) - .5*self.lv) + np.sign(w)
        return w

    def fit(self):
        u = np.random.sample(self.p)
        v = np.random.sample(self.q)
        K = self._K()
        i = 0
        
        while  (np.linalg.norm(u, ord=1) > self.c1) and \
               (np.linalg.norm(v, ord=1) > self.c2):
            u = np.dot(K, v)
            u = self._norm(u)
            u = self._soft_thresh(u, 'u')
            u = self._norm(u)
            v = np.dot(K.T, u)
            v = self._norm(v)
            v = self._soft_thresh(v, 'v')
            v = self._norm(v)
            i += 1
            if i > 100000:
                break
        
        self.K = K
        self.u = u
        self.v = v
        self.i = i

        return self
    
    def predict(self, X_new):
        pass
