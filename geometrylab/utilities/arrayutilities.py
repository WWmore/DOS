#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import sparse

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'


# def sum_repeated_2(array, field):
#     field = field[:]
#     imap = np.zeros(np.amax(field) + 1, dtype=int)
#     index = np.arange(array.shape[0])
#     k, j = np.unique(field, True)
#     imap[k] = np.arange(k.shape[0])
#     result = array[j]
#     field = np.delete(field, j)
#     index = np.delete(index, j)
#     while field.shape[0] > 0:
#         _, j = np.unique(field, True)
#         result[imap[field[j]]] += array[index[j]]
#         field = np.delete(field, j)
#         index = np.delete(index, j)
#     return result


def sum_repeated(array, field):
    N = np.max(field) + 1
    if len(array.shape) == 1:
        j = np.zeros(len(array))
        s = sparse.coo_matrix((array, (field, j)), shape=(N, 1))
        return s.toarray().flatten()
    elif len(array.shape) == 2:
        i = np.tile(field, array.shape[1])
        j = np.repeat(np.arange(array.shape[1]), len(array))
        s = sparse.coo_matrix((array.flatten('F'), (i, j)), shape=(N, array.shape[1]))
        return s.toarray()
    elif len(array.shape) == 3:
        data = np.reshape(array, (array.shape[0], -1))
        i = np.tile(field, data.shape[1])
        j = np.repeat(np.arange(data.shape[1]), len(data))
        s = sparse.coo_matrix((data.flatten('F'), (i, j)), shape=(N, data.shape[1]))
        return np.reshape(s.toarray(), (N, array.shape[1], array.shape[2]))
    else:
        raise NotImplementedError


def repeated_range(array, offset=0):
    if array.shape[0] == 0:
        return np.array([])
    k = np.unique(array)
    imap = np.zeros(np.amax(k) + 1, dtype=int)
    imap[k] = np.arange(offset, offset + k.shape[0])
    rrange = imap[array]
    return rrange


def orthogonal_vectors(array):
    array = format_array(array)
    if len(array.shape) == 1:
        return orthogonal_vector(array)
    ortho = np.zeros(array.shape)
    ortho[:, 0] = - array[:, 1]
    ortho[:, 1] = array[:, 0]
    ortho[np.where((ortho[:, 0] == 0) & (ortho[:, 1] == 0))[0], 1] = 1
    ortho = normalize(ortho)
    return ortho


def orthogonal_vector(array):
    if array[0] == 0 and array[1] == 0:
        if array[2] == 0:
            raise ValueError('zero vector')
        return np.array([1, 0, 0])
    orthogonal = np.array([-array[1], array[0], 0])
    return orthogonal / np.linalg.norm(orthogonal)


def normalize(array, axis=1):
    if len(array.shape) == 1:
        axis = 0
    norm = (np.linalg.norm(array, axis=axis, keepdims=True))
    n = np.divide(array, norm, out=np.zeros(array.shape), where=norm != 0)
    return n


def remap_values(values, target_range=(0, 1)):
    t_int = target_range[1] - target_range[0]
    s_min = np.min(values)
    s_max = np.max(values)
    s_int = s_max - s_min
    s = t_int / s_int
    remapped = (values - s_min) * s + target_range[0]
    return remapped


def format_array(array, shape=None):
    if array is None:
        return None
    if type(shape) is int:
        shape = (shape,)
    if type(array) == float or type(array) == int:
        array = np.array(array)
    elif type(array) == tuple:
        array = list(array)
    if type(array) == list:
        array = np.array(array)
    if shape is None or array.shape == shape:
        return array
    else:
        array_shape = list(array.shape)
        for i in range(len(shape) - len(array.shape)):
            array_shape = [1] + array_shape
        shape = np.array(shape)
        array_shape = np.array(array_shape)
        shape[shape is None] = array_shape[shape is None]
        # repeats = shape - array_shape + 1
        # array = np.tile(array, repeats)
        cmd = 'array['
        for n in shape:
            cmd += ':{},'.format(n)
        cmd += ']'
        array = eval(cmd)
    return array


def bounding_shape(arrays):
    shape = []
    for array in arrays:
        if array is not None:
            if type(array) != np.ndarray:
                array = np.array(array)
            for i in range(len(array.shape) - len(shape)):
                shape = [1] + shape
            for j in range(len(shape)):
                if array.shape[j] > shape[j]:
                    shape[j] = array.shape[j]
    return shape


if __name__ == '__main__':
    import time
    n = 200000
    m = 30
    a = np.random.random((m * n, 2, 2, 2))
    b = np.repeat(np.arange(n), m)
    t0 = time.perf_counter
    d = sum_repeated(a, b)
    print(time.perf_counter() - t0)
