import numpy as np
import pandas as pd


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

def cos_angle(vector1, vector2):
    angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) *
                                        np.linalg.norm(vector2))
    return angle


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


def read_matrices(matrices, labels, num_traits):
    raw_matrices = pd.read_csv(matrices)
    with open(labels) as f:
        monkey_labels = f.read().splitlines()

    def make_symetric(matrix):
        matrix = np.tril(matrix) + np.tril(matrix, k=-1).transpose()
        return matrix

    node_matrices = {monkey_labels[0]: np.array(raw_matrices.ix[0:num_traits-1, ])}
    for i in range(1, len(monkey_labels)):
        new_matrix = np.array(raw_matrices.ix[i*num_traits:(((i+1)*num_traits)-1), :])
        node_matrices[monkey_labels[i]] = make_symetric(new_matrix)

    return node_matrices


def matrix_correlation(matrix1, matrix2):
    tr = matrix1.shape[0]
    correlation = np.corrcoef(matrix1[np.triu_indices(tr,1)], matrix2[np.triu_indices(tr,1)])[1, 0]
    return correlation
