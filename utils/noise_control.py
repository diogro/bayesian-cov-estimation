#! /usr/bin/env python
from numpy import linalg, zeros, dot, diag_indices_from

def a_m(value, vector):
    mat = zeros((len(value), len(value)))
    mat[diag_indices_from(mat)] = value
    mat = dot(dot(vector, mat), vector.transpose())
    return mat

def last_big_eval_index(eig_vals, cutoff):
    """Selects the number of retained eigenvalues based on second
    derivative variance"""
    len_ev = len(eig_vals)
    last_index = None
    for i in range(len_ev):
        if not last_index and eig_vals[i] > cutoff:
            last_index = i

    return last_index

def noise_control(matrix, cutoff = 0.000001):
    """Extends Matrix P and prints it to output_file"""
    eig_vals, eig_vecs = linalg.eigh(matrix)
    last_ev_index = last_big_eval_index(eig_vals, cutoff)
    eig_vals[0:last_ev_index] = eig_vals[last_ev_index]
    p = a_m(eig_vals, eig_vecs)
    return p

