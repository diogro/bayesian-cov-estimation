import numpy as np


def CalcR2(Matrix):
    tr = Matrix.shape[1]
    x, y = np.asarray(np.invert(np.tri(tr, tr, dtype=bool)),
                      dtype=float).nonzero()
    R2Tot = np.mean(Matrix[x, y] * Matrix[x, y])
    return R2Tot


def eigen_var_calc(Matrix):
    eVal, eVec = np.linalg.eigh(Matrix)
    eigenMax = np.zeros(40, float)
    eigenMax[0] = sum(eVal)
    R2Tot = np.var(eVal) / np.var(eigenMax)
    return R2Tot


def icv_calc(Matrix):
    eVal, eVec = np.linalg.eigh(Matrix)
    icv = np.sqrt(np.var(eVal)) / np.mean(eVal)
    return icv


def cos_angle(vector1, vector2):
    angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) *
                                        np.linalg.norm(vector2))
    return angle


def random_skewers(Matrix1, Matrix2):
    n_vector = 1000
    tr = Matrix1.shape[0]
    rand_vec = np.random.multivariate_normal(np.zeros(tr),
                                             np.identity(tr, float),
                                             n_vector).T
    delta_z1 = np.dot(Matrix1, rand_vec)
    delta_z2 = np.dot(Matrix2, rand_vec)
    rs = 0.
    for i in xrange(n_vector):
        rs += cos_angle(delta_z1[:, i], delta_z2[:, i])
    rs = rs / n_vector
    return rs


def matrix_correlation(Matrix1, Matrix2):
    tr = Matrix1.shape[0]
    x, y = np.asarray(np.invert(np.tri(tr, tr, dtype=bool)),
                      dtype=float).nonzero()
    correlation = np.corrcoef(Matrix1[x, y], Matrix2[x, y])[1, 0]
    return correlation
