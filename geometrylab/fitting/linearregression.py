# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'


def linear_regression(points):
    b = np.sum(points, axis=0) / points.shape[0]
    P = np.copy(points)
    for i in range(points.shape[1]):
        P[:,i] -= b[i]
    C = np.einsum('ij,jk->ik', P.T, P)
    eig = np.linalg.eigh(C)
    u = eig[1][:,np.argsort(np.abs(eig[0]))].T
    return  u, b



def fit_sphere(points):
    A = np.zeros((len(points),4))
    A[:,0] = points[:,0]*2
    A[:,1] = points[:,1]*2
    A[:,2] = points[:,2]*2
    A[:,3] = 1
    f = np.zeros((len(points),1))
    f[:,0] = (points[:,0]**2) + (points[:,1]**2) + (points[:,2]**2)
    M = np.dot(A.T,A)
    b = np.dot(A.T,f)
    #C, residules, rank, singval = np.linalg.lstsq(A,f,rcond=-1)
    C = np.linalg.solve(M,b)
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = np.sqrt(t)
    return radius[0], C[0:3].T
















