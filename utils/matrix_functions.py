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


def random_skewers(matrix1, matrix2, num_vectors = 1000):
    traits = matrix1.shape[0]
    rand_vec = np.random.multivariate_normal(np.zeros(traits),
                                             np.identity(traits, float),
                                             num_vectors).T
    delta_z1 = np.dot(matrix1, rand_vec)
    delta_z2 = np.dot(matrix2, rand_vec)

    ndelta_z1 = delta_z1/np.sqrt((delta_z1*delta_z1).sum(0))
    ndelta_z2 = delta_z2/np.sqrt((delta_z2*delta_z2).sum(0))

    return np.mean(np.diag(np.dot(ndelta_z1.T, ndelta_z2))) 

def orandom_skewers(matrix1, matrix2, num_vectors = 1000):
    traits = matrix1.shape[0]
    rand_vec = np.random.multivariate_normal(np.zeros(traits),
                                             np.identity(traits, float),
                                             num_vectors).T
    delta_z1 = np.dot(matrix1, rand_vec)
    delta_z2 = np.dot(matrix2, rand_vec)

    rs = 0.
    for i in xrange(num_vectors):
        rs += cos_angle(delta_z1[:, i], delta_z2[:, i])
    rs = rs / num_vectors

    return rs



def matrix_correlation(Matrix1, Matrix2):
    tr = Matrix1.shape[0]
    x, y = np.asarray(np.invert(np.tri(tr, tr, dtype=bool)),
                      dtype=float).nonzero()
    correlation = np.corrcoef(Matrix1[x, y], Matrix2[x, y])[1, 0]
    return correlation
