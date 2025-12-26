#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'



def clean_quad_boundary(mesh, threshold=1):
    threshold = threshold*0.005
    H = mesh.halfedges
    h1, h2 = mesh.edge_half_edges()
    V1 = mesh.vertices[H[h1,0]]
    V2 = mesh.vertices[H[h2,0]]
    E0 = (V2-V1)/np.linalg.norm(V2-V1, axis=1, keepdims=True)
    V1p = mesh.vertices[H[H[H[H[h1,3],4],3],0]]
    V2p = mesh.vertices[H[H[H[H[h2,3],4],3],0]]
    E1 = (V1p-V1)/np.linalg.norm(V1p-V1, axis=1, keepdims=True)
    E2 = (V2p-V2)/np.linalg.norm(V2p-V2, axis=1, keepdims=True)
    A1 = np.abs(np.einsum('ij,ij->i',E0,E1)) > 1 - threshold
    A2 = np.abs(np.einsum('ij,ij->i',E0,E2)) > 1 - threshold
    C1 = np.logical_or(A1,A2)
    V2p = mesh.vertices[H[H[H[H[H[h1,2],4],2],2],0]]
    V1p = mesh.vertices[H[H[H[H[H[h2,2],4],2],2],0]]
    E1 = (V1p-V1)/np.linalg.norm(V1p-V1, axis=1, keepdims=True)
    E2 = (V2p-V2)/np.linalg.norm(V2p-V2, axis=1, keepdims=True)
    A1 = np.abs(np.einsum('ij,ij->i',E0,E1)) > 1 - threshold
    A2 = np.abs(np.einsum('ij,ij->i',E0,E2)) > 1 - threshold
    C2 = np.logical_or(A1,A2)
    A = np.invert(np.logical_or(C1,C2))
    #----------------------------------
    f1, f2 = mesh.edge_faces()
    L = mesh.face_lengths()
    R = L[f1] + L[f2] - 2
    C3 = R<=5
    boundary = mesh.are_boundary_faces()
    bf1 = boundary[f1]
    bf2 = boundary[f2]
    C4 = np.logical_xor(bf1,bf2)
    B = np.logical_and(C3,C4)
    #-----------------------------------
    v1, v2 = mesh.edge_vertices()
    _,_,R = mesh.vertex_ring_vertices_iterators(sort=True, return_lengths=True)
    C5 = R[v1]+R[v2] > 8
    B = np.logical_and(B,C5)
    A1 = np.logical_and(np.invert(bf1),L[f1]==3)
    A2 = np.logical_and(np.invert(bf2),L[f2]==3)
    C6 = np.logical_or(A1,A2)
    A = np.logical_and(A,C6)
    edges = np.arange(mesh.E)[np.logical_and(A,B)]
    for e in edges:
        mesh.delete_edge(e)
    mesh.clean()

def add_random_noise(mesh, factor=0.1):
    N = mesh.vertex_normals()
    f = np.random.random(mesh.V)
    f -= 0.5
    N[:,0] *= f
    N[:,1] *= f
    N[:,2] *= f
    mesh.vertices += N*factor