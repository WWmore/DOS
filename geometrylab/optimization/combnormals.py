#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import sparse

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.optimization.guidedprojectionbase import GuidedProjectionBase

# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class CombNormals(GuidedProjectionBase):

    _polyline = None

    def __init__(self):
        GuidedProjectionBase.__init__(self)

        weights = {

        'geometric' : 1,

        }

        self.add_weights(weights)

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    @property
    def polyline(self):
        return self._polyline

    @polyline.setter
    def polyline(self, polyline):
        self._polyline = polyline
        self.initialization()

    #--------------------------------------------------------------------------
    #                               Initialization
    #--------------------------------------------------------------------------

    def on_initialize(self):
        tangents = self.polyline.vertex_tangents()
        self.add_value('tangents', np.copy(tangents))
        normals, cos = self.polyline.vertex_bisectors(return_cosines=True)
        self.add_value('normals_0', np.copy(normals))
        self.add_value('normals', np.copy(normals))
        w = np.abs(np.abs(cos) - 1)
        w[w > 0.1] = 1
        if not self.polyline.closed:
            w[0] = w[-1] = 0
        else:
            w[0] = w[-1] = 0
        self.add_value('local_w', np.copy(w))

    def set_dimensions(self):
        self._N = 4*self.polyline.V

    def initialize_unknowns_vector(self):
        normals = self.get_value('normals')
        X = normals.flatten('F')
        X = np.hstack((X, np.ones(self.polyline.V)))
        self._X = X
        self._X0 = np.copy(X)

    def make_errors(self):
        pass

    def post_iteration_update(self):
        pass

    def on_reinitilize(self):
        pass

    #--------------------------------------------------------------------------
    #                                   Utilities
    #--------------------------------------------------------------------------

    def normals(self, untwist=False):
        N = np.reshape(self.X[0:3*self.polyline.V], (self.polyline.V,3), order='F')
        N = utilities.normalize(N)
        vertices = self.polyline.vertices
        if False:#not self.polyline.closed:
            T0 = (vertices[1] - vertices[0])
            N0 = np.cross(np.cross(T0,N[1,:]), T0)
            norm = np.linalg.norm(N0)
            N[0,:] = N0/norm
            T1 = (vertices[-1] - vertices[-2])
            N1 = np.cross(-np.cross(T0,N[-2,:]), T1)
            norm = np.linalg.norm(N1)
            N[-1,:] = N1/norm
        if untwist:
            V0 = N[1:self.polyline.V]
            V1 = N[0:self.polyline.V-1]
            D = np.einsum('ij,ij->i', V0, V1)
            sign = 1
            for i in range(len(D)):
                N[i] *= sign
                if D[i] < -0.2:
                    sign *= -1
            N[-1] = N[-2]
        return N

    def plot_normals(self):
        from geometrylab import vtkplot
        N = vtkplot.Vectors(self.normals(), anchor=self.polyline.vertices)
        w = self.get_value('local_w')
        w = np.column_stack((w,w,w))
        N0 = vtkplot.Vectors(self.get_value('normals')*w,
                             anchor=self.polyline.vertices,
                             color='gray_50')
        #N00 = vtkplot.Vectors(self.get_value('normals_0'), anchor=self.polyline.vertices, color='b')
        vtkplot.view([N,N0])
    # -------------------------------------------------------------------------
    #                          Geometric Constraints
    # -------------------------------------------------------------------------

    def unit_constraints(self):
        w = 2
        V = self.polyline.V
        v = np.arange(V)
        i = np.arange(V)
        i = np.hstack((i, i, i))
        j = np.hstack((v, V+v, 2*V+v))
        X = self.X
        data = 2 * np.hstack((X[v], X[V+v], X[2*V+v])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(V, self.N))
        r = ((X[v]**2 + X[V+v]**2 + X[2*V+v]**2) + 1) * w
        self.add_iterative_constraint(H, r, 'unit')

    def tangent_constraints(self):
        w = 10 #- self.get_value('local_w')
        V = self.polyline.V
        t = self.get_value('tangents')
        v = np.arange(V)
        i = np.arange(V)
        i = np.hstack((i, i, i))
        j = np.hstack((v, V+v, 2*V+v))
        data = np.hstack((t[:,0]*w, t[:,1]*w, t[:,2]*w))
        H = sparse.coo_matrix((data,(i,j)), shape=(V, self.N))
        r = np.zeros(V)
        self.add_constant_constraint(H, r, 'tangent')

    def fixed_constraints(self):
        w = self.get_value('local_w') * 1
        V = self.polyline.V
        n = self.get_value('normals')
        v = np.arange(V)
        one = w * np.ones(V)
        i = np.hstack((v, V+v, 2*V+v))
        j = np.hstack((v, V+v, 2*V+v))
        data = np.hstack((one, one, one))
        r = n.flatten('F') * np.hstack((w,w,w))
        H = sparse.coo_matrix((data,(i,j)), shape=(3*V, self.N))
        self.add_constant_constraint(H, r, 'fixed')

    def parallel_constraints(self):
        w = self.get_value('local_w') * 1
        V = self.polyline.V
        n = self.get_value('normals')
        v = np.arange(V)
        i = np.hstack((v, V+v, 2*V+v))
        j = np.hstack((v, V+v, 2*V+v))
        d = self.X[3*V:4*V] * w
        data = np.hstack((d, d, d))
        r = n.flatten('F') * np.hstack((w,w,w))
        H = sparse.coo_matrix((data,(i,j)), shape=(3*V, self.N))
        self.add_iterative_constraint(H, r, 'parallel')

    def normals_fairness(self):
        w = 0.1
        v0 = np.arange(self.polyline.V)
        v2 = np.roll(v0,1)
        v1 = np.roll(v0,-1)
        if not self.polyline.closed:
            v0 = v0[1:-1]
            v2 = v2[1:-1]
            v1 = v1[1:-1]
        W = len(v0)
        V = self.polyline.V
        one  = np.ones(W)
        i = np.arange(W)
        i = np.hstack((i,i,i))
        j = np.hstack((v0,v1,v2))
        data = np.hstack((-2*one, one, one)) * w
        i = np.hstack((i, W+i, 2*W+i))
        j = np.hstack((j, V+j, 2*V+j))
        data = np.hstack((data, data, data))
        K = sparse.coo_matrix((data,(i,j)), shape=(3*W, self.N))
        s = np.zeros(3*W)
        self.add_constant_fairness(K, s, 'fairness')

    # -------------------------------------------------------------------------
    #                                 Build
    # -------------------------------------------------------------------------

    def build_iterative_constraints(self):
        self.unit_constraints()
        self.tangent_constraints()
        #self.fixed_constraints()
        self.parallel_constraints()

    def build_constant_fairness(self):
        self.normals_fairness()

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def comb_normals(polyline, untwist=False):
    optimizer = CombNormals()
    optimizer.verbose = False
    optimizer.polyline = polyline
    optimizer.iterations = 5
    optimizer.epsilon = 0.0001
    optimizer.optimize()
    #optimizer.plot_normals()
    return optimizer.normals(untwist)
