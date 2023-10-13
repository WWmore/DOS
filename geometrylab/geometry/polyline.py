#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import sparse

from scipy import spatial

from scipy.sparse import linalg

# -----------------------------------------------------------------------------

from geometrylab import utilities

# -----------------------------------------------------------------------------
'''Hui add: cell_array(diagonal polylines) '''
'''polyine.py: The Polyline data structure class'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                   POLYLINE
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Polyline(object):

    def __init__(self, vertices, closed=False):

        self.name = 'polyline'

        if type(vertices) == str:
            self.read_obj_file(vertices)
        else:
            self.vertices = vertices

        self.closed = closed

        self.corner_tolerance = None

        self._kdtree = None

        self._cell_array = None # Hui
    #--------------------------------------------------------------------------
    #                               Attributes
    #--------------------------------------------------------------------------

    @property
    def type(self):
        return 'Polyline'

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def E(self):
        if self.closed:
            return self.V
        else:
            return self.V - 1

    @property
    def vertices(self):
        return self._vertices

    @vertices.setter
    def vertices(self, vertices):
        vertices = np.array(vertices, 'f')
        if len(vertices.shape) != 2 or vertices.shape[1] != 3:
            raise ValueError('wrong size!')
        try:
            self._vertices[:,:] = vertices[:,:]
        except:
            self._vertices = vertices
        self._kdtree = None

    @property
    def corner_tolerance(self):
        return self._corner_tolerance

    @corner_tolerance.setter
    def corner_tolerance(self, corner_tolerance):
        if corner_tolerance is None:
            self._corner_tolerance = None
        else:
            self._corner_tolerance = float(corner_tolerance)

    @property
    def closed(self):
        return self._closed

    @closed.setter
    def closed(self, bool):
        self._closed = bool

    @property
    def cell_array(self): # Hui add: not the one used in polylinesource.py
        return self._cell_array

    @cell_array.setter
    def cell_array(self, cells): # Hui add: write by other functions
        self._cell_array = cells
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------

    def __str__(self):
        if self.closed:
            out = 'Closed '
        else:
            out = 'Open '
        out += ('polyline: |V| = {}').format(self.V)

        return out

    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------

    def read_obj_file(self, file_name):
        file_name = str(file_name)
        self.name = file_name.split('.')[0]
        obj_file = open(file_name, encoding='utf-8')
        vertices_list = []
        for l in obj_file:
            splited_line = l.split(' ')
            if splited_line[0] == 'v':
                split_x = splited_line[1].split('\n')
                x = float(split_x[0])
                split_y = splited_line[2].split('\n')
                y = float(split_y[0])
                split_z = splited_line[3].split('\n')
                try:
                    z = float(split_z[0])
                except ValueError:
                    print('WARNING: disable line wrap when saving .obj')
                vertices_list.append([x, y ,z])
        self.vertices = np.array(vertices_list)


    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------

    def vertex_tangents(self):
        i = np.arange(self.V)
        v1 = np.roll(i,1)
        v2 = np.roll(i,-1)
        if not self.closed:
            v1[0] = 0
            v2[-1] = self.V - 1
        N = self.vertices[v2,:] - self.vertices[v1,:]
        N = N / np.linalg.norm(N, axis=1, keepdims=True)
        return N

    def vertex_bisectors(self, return_cosines=False):
        v0 = np.arange(self.V)
        v2 = np.roll(v0,1)
        v1 = np.roll(v0,-1)
        V1 = self.vertices[v1,:] - self.vertices[v0,:]
        V2 = self.vertices[v0,:] - self.vertices[v2,:]
        V1 = utilities.normalize(V1)
        V2 = utilities.normalize(V2)
        B = V1 - V2
        B = utilities.normalize(B)
        if not self.closed:
            N0 = np.cross(np.cross(V1[0],B[1,:]), V1[0])
            norm = np.linalg.norm(N0)
            B[0,:] = N0/norm
            N1 = np.cross(np.cross(V2[-1],B[-2,:]), V2[-1])
            norm = np.linalg.norm(N0)
            B[-1,:] = N1/norm
        if not return_cosines:
            return B
        else:
            A = np.einsum('ij,ij->i', V1, -V2)
            if not self.closed:
                A[0] = A[-1] = -1
            return B, A

    def corners(self):
        if self.corner_tolerance is None:
            return []
        v1 = np.arange(self.V)
        v0 = np.roll(v1,1)
        v2 = np.roll(v1,-1)
        if not self.closed:
            v0[0] = self.V - 1
            v2[-1] = 0
        V1 = self.vertices[v2,:] - self.vertices[v1,:]
        V1 = utilities.normalize(V1)
        V2 = self.vertices[v1,:] - self.vertices[v0,:]
        V2 = utilities.normalize(V2)
        C = np.einsum('ij,ij->i', V1, V2)
        corners = v1[np.where(C[:] < self.corner_tolerance)[0]].tolist()
        return corners

    def refine(self, steps=5):
        if steps <= 0:
            return
        V = self.V
        corners = self.corners()
        if self.closed:
            N = steps*V
        else:
            N = (steps)*(V-1)
        v = np.arange(V)
        s = np.arange(V+N)
        d  = np.repeat(3, V+N)
        d1 = np.repeat(-2, V+N)
        d2 = np.repeat(0.5, V+N)
        d[(steps+1)*v] = 1
        d1[(steps+1)*v] = 0
        d2[(steps+1)*v] = 0
        i = np.hstack((s,s,s,s,s))
        j = np.hstack((s, np.roll(s,1), np.roll(s,-1),
                       np.roll(s,2), np.roll(s,-2)))
        data = np.hstack((d,d1,d1,d2,d2))
        if not self.closed:
            data[1] = data[V+N-2] = 2.5
            data[V+N+1] = data[3*(V+N)-2] = -1
            data[3*(V+N)+1] = 0
            data[5*(V+N)-2] = 0
            if corners is not None and len(corners) > 0:
                corners = np.array(corners)
                corners = np.delete(corners, np.where(corners==0)[0])
        if corners is not None and len(corners) > 0:
            c = (steps+1)*np.array(corners)
            data[np.roll(s,1)[c]] = 2.5
            data[np.roll(s,-1)[c]] = 2.5
            data[V+N+np.roll(s,-1)[c]] = -1
            data[2*(V+N)+np.roll(s,1)[c]] = -1
            data[3*(V+N)+np.roll(s,-1)[c]] = 0
            data[4*(V+N)+np.roll(s,1)[c]] = 0
        M = sparse.coo_matrix((data,(i,j)), shape=(V+N,V+N))
        P = np.zeros((N+V,3))
        P[(steps+1)*v] = self.vertices
        M = sparse.csc_matrix(M)
        X = linalg.spsolve(M,P[:,0], use_umfpack=False)
        Y = linalg.spsolve(M,P[:,1], use_umfpack=False)
        Z = linalg.spsolve(M,P[:,2], use_umfpack=False)
        P = np.vstack((X,Y,Z)).T
        self.vertices = P

#------------------------------------------------------------------------------
#                                  CLOSENESS
#------------------------------------------------------------------------------

    def make_kdtree(self):
        KDtree = spatial.cKDTree(self.vertices)
        self._kdtree = KDtree

    def closest_vertices(self, points, make_tree=False):
        if self._kdtree is None:
            self.make_kdtree()
        elif make_tree:
            self.make_kdtree()
        closest = self._kdtree.query(points)[1]
        return closest

#------------------------------------------------------------------------------
#                                VISUALIZATION
#------------------------------------------------------------------------------

    def cell_array(self, offset=0):
        vi = np.arange(self.V) + offset
        vj = np.roll(vi,-1) + offset
        c = np.repeat(2,self.vertices.shape[0])
        cells = np.vstack((c,vi,vj)).T
        if not self.closed:
            cells = np.delete(cells, self.V-1, axis=0)
        cells = np.ravel(cells)
        return cells
