import numpy as np
import gc

def rand_corr(n, ke):
    """Return a n x n random correlation matrix"""
    c = np.zeros((n,n))
    b = np.tri(n, n)

    c[1:n, 0] = -1 + 2*np.round(np.random.rand(1, n-1)*10**8)/10**8
    b[1:n, 0] = c[1:n, 0]

    for i in xrange(1, n):
        b[i, 1:i+1] = np.sqrt(1 - c[i, 0]**2)

    for i in xrange(2, n):
        for j in xrange(1, i):
            b1 = np.dot(b[j, 0:j], b[i, 0:j].T)
            b2 = np.dot(b[j, j], b[i, j])
            z = b1 + b2
            y = b1 - b2
            if b2 < ke:
                c[i, j] = b1
                cosinv = 0
            else:
                c[i, j] = y + (z - y)*np.round(np.random.rand()*10**8)/10**8

            cosinv = (c[i, j] - b1)/b2

            if np.isfinite(cosinv):
                if cosinv > 1:
                    b[i, j] = b[i, j]
                    b[i, j+1:n+1] = 0
                elif cosinv < -1:
                    b[i, j] = -b[i, j]
                    b[i, j+1:n+1] = 0
                else:
                    b[i, j] = b[i, j]*cosinv
                    sinTheta = np.sqrt(1 - cosinv**2)
                    for k in xrange(j+1, n):
                        b[i, k] = b[i, k]*sinTheta

    c = c + c.T + np.eye(n)
    perm = np.random.permutation(n)
    c = (c[perm,])[:,perm]

    return c

def triang_decomp(c):
    """Return a lower triangular matrix B that B * B.T = C"""
    n = c.shape[0]
    b = np.tri(n, n)

    b[1:n, 0] = c[1:n, 0]

    for i in xrange(1, n):
        b[i, 1:i+1] = np.sqrt(1 - c[i, 0]**2)

    for i in xrange(2, n):
        for j in xrange(1, i):
            b1 = np.dot(b[j, 0:j], b[i, 0:j].T)
            b2 = np.dot(b[j, j], b[i, j])
            cosinv = (c[i, j] - b1)/b2

            if np.isfinite(cosinv):
                if cosinv > 1:
                    b[i, j] = b[i, j]
                    b[i, j+1:n+1] = 0
                elif cosinv < -1:
                    b[i, j] = -b[i, j]
                    b[i, j+1:n+1] = 0
                else:
                    b[i, j] = b[i, j]*cosinv
                    sinTheta = np.sqrt(1 - cosinv**2)
                    for k in xrange(j+1, n):
                        b[i, k] = b[i, k]*sinTheta

    return b

def calc_params(b):
    """ Given a lower trianguler B matrix returned by triang_decomp, return angle parameters"""
    n = b.shape[0]
    p = np.zeros(n*n).reshape((n, n))

    for j in xrange(n-1):
        for i in xrange(j, n):
            p[i, j] = np.arccos(b[i, j]/np.exp(np.sum(np.log(np.sin(p[i, 0:j])))))

    p[np.isnan(p)] = 0

    return p

def triang_from_params(p):
    """Given a p param matrix, calculate the lower triangular B matrix"""
    n = p.shape[0]
    b = np.zeros(n*n).reshape((n,n))

    for j in xrange(n-1):
        for i in xrange(n):
            b[i, j] = np.cos(p[i,j])*np.product(np.sin(p[i, 0:j]))

    b[n-1,n-1] = np.product(np.sin(p[n-1, 0:n-1]))

    return b

def calc_path(ms, bp, i, j, s=100, pcs=[0], calc_rs=False, calc_ev=False, return_matrices=False):
    m0 = ms[i]
    m0_initial = ms[i]
    r2 = np.zeros(s)
    flex = np.zeros(s)
    isoc = []
    evs = []
    rs = []

    if return_matrices:
        ms_int = [m0]

    iso = np.ones(m0.shape[0])/np.sqrt(m0.shape[0])
    diff = (bp[j][1] - bp[i][1])/100
    p = bp[i][1]

    for i in xrange(100):
        if calc_ev:
            evs.append(np.linalg.eigh(m0)[0])
        r2[i] = calc_r2(m0)
        flex[i] = flexibility(m0)
        for j in pcs:
            isoc.append(np.abs(np.dot(np.linalg.eig(m0)[1][:,j], iso)))

        if calc_rs:
            rs.append(random_skewers(m0_initial, m0))

        p += diff
        new_b = triang_from_params(p)
        del m0
        m0 = np.dot(new_b, new_b.T)

        if return_matrices:
            ms_int.append(m0)

    gc.collect()
    if calc_rs:
        return r2, flex, isoc, rs

    if not calc_ev:
        if return_matrices:
            return r2, flex, isoc, ms_int
    else:
        return r2, flex, isoc, evs

def random_skewers(matrix1, matrix2, num_vectors=1000):
    traits = matrix1.shape[0]
    rand_vec = np.random.multivariate_normal(np.zeros(traits),
                                             np.identity(traits, float),
                                             num_vectors).T
    delta_z1 = np.dot(matrix1, rand_vec)
    delta_z2 = np.dot(matrix2, rand_vec)

    ndelta_z1 = delta_z1/np.sqrt((delta_z1*delta_z1).sum(0))
    ndelta_z2 = delta_z2/np.sqrt((delta_z2*delta_z2).sum(0))

    return np.mean(np.diag(np.dot(ndelta_z1.T, ndelta_z2)))

def flexibility(matrix1, num_vectors=1000):
    traits = matrix1.shape[0]
    rand_vec = np.random.multivariate_normal(np.zeros(traits),
                                             np.identity(traits, float),
                                             num_vectors).T

    rand_vec = rand_vec/np.sqrt((rand_vec*rand_vec).sum(0))

    delta_z1 = np.dot(matrix1, rand_vec)

    ndelta_z1 = delta_z1/np.sqrt((delta_z1*delta_z1).sum(0))

    return np.mean(np.diag(np.dot(ndelta_z1.T, rand_vec)))

def calc_r2(m):
   tr = m.shape[1]
   x, y = np.asarray(np.invert(np.tri(tr, tr, dtype=bool)), dtype=float).nonzero()
   r2_tot = np.mean(m[x, y] * m[x, y])
   return r2_tot

