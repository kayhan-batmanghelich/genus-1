from __future__ import division
import numpy as np
"""
FUNCTIONS NOT INCORPORATED FROM fixedformVB.m:

function   outIdx = item2Feature(inIdx, itemToFeatureMap)
        outIdx = [] ;
        for idx = inIdx(:)'
            outIdx = union(outIdx,itemToFeatureMap(idx).idx) ;
        end
end

"""
def elem_sympoly(lamb, k):
    N = len(lamb)
    E = np.zeros((k+1, N+1))
    for l in np.arange(k)+1:
        for n in np.arange(N)+1:
            E[l, n] = E[l, n-1]+(lamb[n-1]*E[l-1,n-1])
    return E
            
def sample_k(lamb, k):
    E = elem_sympoly(lamb, k)
    i = len(lamb)
    rem = k
    S = np.zeros((k, 1))
    while rem>0:
        if i == rem:
            marg = 1
        else:
            marg = lamb(i)*E[rem, i]/E[rem, i+1]
        if np.random.random() < marg:
            S[rem] = i
            rem -= 1
        i -= 1
    return S

def sample_dual_dpp(B, C, k):
    if k:
        D = C['D']/(1+C['D'])
        v = np.where(np.random.random(size=(len(D), 1)) <= D)
    else:
        v = sample_k(C['D'], k)
    k = len(v)
    V = C['V'][:, v]  
    V = V*1./(np.sqrt(C['D'][v].T))  
    Y = np.zeros((k, 1))
    i_iter = np.arange(k)[::-1]
    for i in i_iter:
        P = np.sum(np.dot(B, V)**2, axis=1)
        P[Y[np.arange(k, i+1)[::-1]]] = 0
        P = P/np.sum(P)
        Y[i] = np.where(np.random.random() <= np.cumsum(P))[0][0] # not sure about this
        S = np.dot(B[Y[i], :], V)
        j = np.where(S)[0]
        Vj = V[:, j]
        Sj = S[j]
        # equiv matlab code:
        # V = V(:,[1:j-1 j+1:end]);
        V = V[:, (np.arange(j-1), np.arange(j+1,i_iter[:, i][-1]))] # not sure about this
        S = S[:, (np.arange(j-1), np.arange(j+1,i_iter[:, i][-1]))]
        V = V-(Vj*(S/Sj))
        for a in np.arange(i-1):
            for b in np.arange(a-1):
                V[:, a] = V[:, a]-(np.dot(V[:, a].T, C['M']).dot(V[:, b])).dot(V[:, b])
            V[:, a] = V[:, a]/np.sqrt(np.dot(V[:, a].T, C['M']).dot(V[:, a]))
        if np.any(np.isnan(V)):
            Y = Y[Y>0]
            print("There are NaNs in V")
    Y = np.sort(Y)
    if len(Y) != len(np.unique(Y)):
        print("Redundancy in the set")
        Y = np.unique(Y)
    return Y

def decompose_kernel(M):
    L = {'M':np.copy(M)}
    V, D = np.eig(M)
    L['V'] = np.real(V)
    L['D'] = np.real(np.diag(D))
    return L

def compute_Khelper(Phi, eta):
    N = np.shape(Phi)[0]
    L = np.dot(np.diag(np.exp(eta/2)), Phi)
    L = np.dot(L, L.T)
    return np.linalg.solve((L + np.eye(N)), L)

def sample_from_DPP(eta, Phi, opt):
    status = 0
    try:
        q = np.exp(eta)
        C = decompose_kernel(Phi.T.dot(np.diag(q**2)).dot(Phi))
    except:
        status = 1
    num_items = np.shape(Phi)[0]
    S = np.empty(opt['num_batch_samp']) # not sure about this
    S_matrix = np.zeros(num_items, opt['num_batch_samp'])   
    if np.any(C['D']<=0):
        print("At least one of the eigenvalues is negative, sign of numerical instability")
        idx = np.logical_or(C['D']<=0, np.imag(C['D'])!=0)
        C['D'][idx] = 0
        C['V'][:, idx] = 0     
    for cnt in np.arange(opt['num_batch_samp']):
        if opt['use_KDpp']:
            S[cnt] = sample_dual_dpp(np.dot(np.diag(q), Phi), C, opt['exp_card'])
        else:
            S[cnt] = sample_dual_dpp(np.dot(np.diag(q), Phi), C)
        S_matrix[S[cnt], cnt] = 1
    return S_matrix, status

def compute_inclustion_proba(y, X, S_matrix, opt):
# not implemented yet
    pass
    
def fixed_formVB(y, X, Phi, opt, item_to_fmap=None):
    if not item_to_fmap:
        item_to_fmap = []
    max_iter = opt['max_iter']
    num_items = np.shape(Phi)[0]
    if opt['num_draw_before_inv']:
        num_draw_before_inv = opt['num_draw_before_inv']
    else:
        num_draw_before_inv = 1
    if opt['use_intercept']:
        use_intercept = opt['use_intercept']
    else:
        use_intercept = False
    if max_iter < num_draw_before_inv:
        raise Exception("Max iterations should be larger than 3*num_items")
    if opt['initial_fcn']:
        reg_weight = 0
        eta, C = opt['initial_Fcn']
    else:
        if opt['normalize_Phi']:
            # matlab code:
            # q0 = sqrt( sum(Phi.*conj(Phi) ,2))
            q0 = np.sqrt(np.sum(Phi*Phi.conj(), axis=1)) # not sure about this, might be a dot product
            Phi = Phi*(1./q0)
            Phi[q0==0, :] = 0
        else:
            q0 = np.ones((num_items, 1))
        L = np.dot(np.diaq(q0), Phi.dot(Phi.T)).dot(np.diag(q0))
        reg_weight = 0.
        if not opt['stop_eta0']:
            opt['stop_eta0'] = np.inf
        if not opt['plot_eta']:
            opt['plot_eta'] = True
        # matlab code:
        # eta = [2*log(q0) ; 0];
        eta = np.array([2*np.log(q0), 0])[:, None] # not sure about this
        # matlab code:
        # K = computeKHelper(Phi, eta(1:end-1));
        K = compute_Khelper(Phi, eta[:-1]) # not sure about this
        if opt['use_marginal_4C']:
            K_last_col_p1 = np.array(np.diag(K).T.tolist()+[1])
            C = np.vstack((np.hstack((K, np.diag(K))), K_last_col_p1))
        else:
            C_top = np.hstack((np.diag(np.diag(k)), np.zeros((num_items, 1))))
            C_bottom = np.array(np.zeros((1, num_items)).tolist()+[1])
            C = np.vstack((C_top, C_bottom))
    g = C*eta
    C_bar = np.zeros((num_items+1, num_items+1))
    g_bar = np.zeros((num_items+1, 1))
    w = 1/np.sqrt(max_iter)
    # not done yet
