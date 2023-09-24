#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'


def sum_repeated(array, field):
    field = field[:]
    imap = np.zeros(np.amax(field)+1, dtype=np.int)
    index = np.arange(array.shape[0])
    k, j = np.unique(field, True)
    imap[k] = np.arange(k.shape[0])
    result = array[j]
    field = np.delete(field,j)
    index = np.delete(index,j)
    while field.shape[0] > 0:
        _, j = np.unique(field, True)
        result[imap[field[j]]] += array[index[j]]
        field = np.delete(field,j)
        index = np.delete(index,j)
    return result

def repeated_range(array, offset=0):
    if array.shape[0] == 0:
        return np.array([])
    k = np.unique(array)
    imap = np.zeros(np.amax(k) + 1, dtype=np.int)
    imap[k] = np.arange(offset, offset + k.shape[0])
    rrange = imap[array]
    return rrange

def orthogonal_vectors(array):
    O = np.zeros(array.shape)
    O[:,0] = - array[:,1]
    O[:,1] = array[:,0]
    O[np.where((O[:,0] == 0) & (O[:,1] == 0))[0],1] = 1
    O = O / np.linalg.norm(O, axis=1, keepdims=True)
    return O

def orthogonal_vector(array):
    if array[0] == 0 and array[1] == 0:
            if array[2] == 0:
                raise ValueError('zero vector')
            return np.array([1,0,0])
    O = np.array([-array[1],array[0],0])
    return O / np.linalg.norm(O)

def normalize(array, axis=1):
    eps = 1e-10
    array = array / (np.linalg.norm(array, axis=axis, keepdims=True) + eps)
    return array

def remap(source, target_range=(0,1)):
    t_int = target_range[1] - target_range[0]
    s_min = np.min(source)
    s_max = np.max(source)
    s_int = s_max - s_min
    s = t_int/s_int
    remapped = (source - s_min) * s + target_range[0]
    return remapped


if __name__ == '__main__':
    a = np.array([0,8,2,1,1,3,4,6,6,8,8,8])
    b = repeated_range(a)
    print(b)



