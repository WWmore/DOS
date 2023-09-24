# -*- coding: utf-8 -*-
from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy.sparse import coo_matrix

from pypardiso import spsolve



# -----------------------------------------------------------------------------

from geometrylab import utilities

# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------

def quadratic_fitting(mesh, double_ring=False, interpolate_vertex=False):
    e3 = mesh.vertex_normals()
    e1 = utilities.orthogonal_vectors(e3)
    e2 = np.cross(e3,e1)
    if double_ring:
        vi, vj, = mesh.vertex_double_ring_vertices_iterators()
    else:
        vi, vj, = mesh.vertex_ring_vertices_iterators()
    N = vj.shape[0]
    Pj = mesh.vertices[vj] - mesh.vertices[vi]
    x1 = np.einsum('ij,ij -> i', Pj, e1[vi])
    x2 = np.einsum('ij,ij -> i', Pj, e2[vi])
    x3 = np.einsum('ij,ij -> i', Pj, e3[vi])
    i = np.arange(N)
    if interpolate_vertex:
        i = np.hstack((i,i,i,i,i))
        j = 5*vi
        j = np.hstack((j,j+1,j+2,j+3,j+4))
        data = np.hstack((x1, x2, x1**2, x1*x2, x2**2))
        D = coo_matrix((data, (i, j)), shape=(N, 5*mesh.V))
        D = D.tocsr()
        A = (D.T).dot(D)
        b = (D.T).dot(x3)
        a = spsolve(A, b)
        a = np.reshape(a, (mesh.V, 5))
        a = np.column_stack((np.zeros(mesh.V), a))
    else:
        i = np.hstack((i,i,i,i,i,i))
        j = 6*vi
        j = np.hstack((j,j+1,j+2,j+3,j+4,j+5))
        one = np.ones(x1.shape[0])
        data = np.hstack((one, x1, x2, x1**2, x1*x2, x2**2))
        D = coo_matrix((data, (i, j)), shape=(N, 6*mesh.V))
        D = D.tocsr()
        A = (D.T).dot(D)
        b = (D.T).dot(x3)
        a = spsolve(A, b)
        a = np.reshape(a, (mesh.V, 6))
    return a

def principal_curvatures(mesh, double_ring=False, interpolate_vertex=False):
    e3 = mesh.vertex_normals()
    e1 = utilities.orthogonal_vectors(e3)
    e2 = np.cross(e3, e1)
    a = quadratic_fitting(mesh, double_ring, interpolate_vertex)
    g = (a[:,1]**2 + a[:,2]**2) + 1
    W = np.array([[2*g*a[:,3], g*a[:,4]],[g*a[:,4], 2*g*a[:,5]]]).T
    eig = np.linalg.eigh(W)
    k1 = eig[0][:,0]
    k2 = eig[0][:,1]
    V1 = eig[1][:,:,0]
    V2 = eig[1][:,:,1]
    V1 = (e1.T * V1[:,0]).T + (e2.T * V1[:,1]).T
    V2 = (e1.T * V2[:,0]).T + (e2.T * V2[:,1]).T
    return (k1, k2, V1, V2)

'''
def principal_curvatures2(mesh, double_ring=False, interpolate_vertex=True):
    e3 = mesh.vertex_normals()
    e1 = utilities.orthogonal_vectors(e3)
    e2 = np.cross(e3,e1)
    if double_ring:
        vi, vj, = mesh.vertex_double_ring_vertices_iterators()
    else:
        vi, vj, = mesh.vertex_ring_vertices_iterators()
    N = vj.shape[0]
    Pj = mesh.vertices[vj] - mesh.vertices[vi]
    x1 = np.einsum('ij,ij -> i', Pj, e1[vi])
    x2 = np.einsum('ij,ij -> i', Pj, e2[vi])
    x3 = np.einsum('ij,ij -> i', Pj, e3[vi])
    i = np.arange(N)
    if interpolate_vertex:
        i = np.hstack((i,i,i,i,i))
        j = 5*vi
        j = np.hstack((j,j+1,j+2,j+3,j+4))
        data = np.hstack((x1, x2, x1**2, x1*x2, x2**2))
        D = coo_matrix((data, (i, j)), shape=(N, 5*mesh.V))
        D = D.tocsr()
        A = (D.T).dot(D)
        b = (D.T).dot(x3)
        a = spsolve(A, b)
        a1 = a[::5]
        a2 = a[1::5]
        a3 = a[2::5]
        a4 = a[3::5]
        a5 = a[4::5]
    else:
        i = np.hstack((i,i,i,i,i,i))
        j = 6*vi
        j = np.hstack((j,j+1,j+2,j+3,j+4,j+5))
        one = np.ones(x1.shape[0])
        data = np.hstack((one, x1, x2, x1**2, x1*x2, x2**2))
        D = coo_matrix((data, (i, j)), shape=(N, 6*mesh.V))
        D = D.tocsr()
        A = (D.T).dot(D)
        b = (D.T).dot(x3)
        a = spsolve(A, b)
        a1 = a[1::6]
        a2 = a[2::6]
        a3 = a[3::6]
        a4 = a[4::6]
        a5 = a[5::6]
    g = (a1**2 + a2**2) + 1
    W = np.array([[2*g*a3, g*a4],[g*a4, 2*g*a5]]).T
    eig = np.linalg.eigh(W)
    k1 = eig[0][:,0]
    k2 = eig[0][:,1]
    V1 = eig[1][:,:,0]
    V2 = eig[1][:,:,1]
    V1 = (e1.T * V1[:,0]).T + (e2.T * V1[:,1]).T
    V2 = (e1.T * V2[:,0]).T + (e2.T * V2[:,1]).T
    return (k1, k2, V1, V2)
'''