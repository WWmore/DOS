#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.optimization.combnormals import comb_normals

# -----------------------------------------------------------------------------

'''facesource.py: The mesh faces plot source class'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                               Mesh Vector Glyph
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class Arrows(object):

    def __init__(self, vectors, anchor, **kwargs):

        self.vectors = vectors

        self._anchor = anchor

        self._shaft_radius = kwargs.get('shaft_radius', None)

        self._sides = kwargs.get('sides', 20)

        self._tip_length = kwargs.get('tip_length', None)

        self._tip_radius = kwargs.get('tip_radius', None)

        self._scale_factor = kwargs.get('scale_factor', 1)

        self._position = kwargs.get('position', 'head')

        self._offset = kwargs.get('offset', 0)

        self._mask = None

        self._make_faces()

        self._make_vertices()

    @property
    def V(self):
        return self.vectors.shape[0]

    @property
    def F(self):
        return self._F

    @property
    def G(self):
        return self._G

    @property
    def lengths(self):
        return self._lengths

    @property
    def vectors(self):
        return self._vectors[self.mask]

    @vectors.setter
    def vectors(self, vectors):
        try:
            vectors.shape
        except AttributeError:
            vectors = np.array(vectors)
        if len(vectors.shape) == 1:
            vectors = np.array([vectors])
        self._vectors = vectors

    @property
    def anchor(self):
        return self._anchor[self._mask]

    @anchor.setter
    def anchor(self, anchor):
        self._anchor = anchor

    @property
    def sides(self):
        return self._sides

    @sides.setter
    def sides(self, sides):
        if self._sides != sides:
            self._sides = sides
            self.make_faces()

    @property
    def mask(self):
        if self._mask is None:
            norm = np.linalg.norm(self._vectors, axis=1)
            self._mask = norm != 0
        return self._mask

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        update_faces = False
        self._mask = None
        self.vectors = kwargs.get('vectors', self.vectors)
        self._anchor = kwargs.get('anchor', self._anchor)
        self._shaft_radius = kwargs.get('shaft_radius', self._shaft_radius)
        if self._sides != kwargs.get('sides', self._sides):
            self._sides = kwargs.get('sides', self._sides)
            update_faces = True
        self._tip_length = kwargs.get('tip_length', self._tip_length)
        self._scale_factor = kwargs.get('scale_factor', self._scale_factor)
        self._position = kwargs.get('position', self._position)
        self._offset = kwargs.get('offset', self._offset)
        if update_faces:
            self._make_faces()
        self._make_vertices()

    def _make_faces(self):
        try:
            self._vectors.shape
        except:
            self._vectors = np.array(self._vectors)
        if len(self._vectors.shape) == 1:
            self._vectors = np.array([self._vectors])
        try:
                self._anchor.shape
        except:
            self._anchor = np.array(self._anchor)
        if len(self._anchor.shape) == 1:
            self._anchor = np.array([self._anchor])
        Fx = self._sides; Fy = 7
        N = self.vectors.shape[0]
        F = N * Fx * Fy + 1
        M = np.arange(Fx*Fy, dtype=int)
        M = np.reshape(M,(Fy,Fx))
        M = np.insert(M,Fx,M[:,0],axis=1)
        Q = np.zeros((Fy-1, Fx,4) ,dtype=int)
        Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
        Q[:,:,2] = M[:M.shape[0]-1,1:]
        Q[:,:,1] = M[1:,1:]
        Q[:,:,0] = M[1:,:M.shape[1]-1]
        T = np.zeros((1, Fx, 3) ,dtype=int)
        T[0,:,0] = M[-1,1:]
        T[0,:,1] = M[-1,:M.shape[1]-1]
        T[0,:,2] = Fx*Fy
        O = np.arange(0,N*(Fx*Fy+1),Fx*Fy+1)
        O = np.reshape(O,(N,1,1))
        Q = Q.reshape(((Fx)*(Fy-1),4),order='C')
        T = T.reshape((Fx,3),order='C')
        t1 = np.repeat(1,N*Fx*(Fy-1))
        t2 = np.repeat(0,N*Fx)
        types = np.hstack((t1,t2))
        Qd = np.zeros((N,Q.shape[0],Q.shape[1]),dtype=int)
        Qd[:] = Q
        Qd += O
        Td = np.zeros((N,T.shape[0],T.shape[1]),dtype=int)
        Td[:] = T
        Td += O
        Qd = np.insert(Qd,0,4,axis=2)
        Td = np.insert(Td,0,3,axis=2)
        Qd = np.reshape(Qd,(Qd.shape[0]*Qd.shape[1]*Qd.shape[2]))
        Td = np.reshape(Td,(Td.shape[0]*Td.shape[1]*Td.shape[2]))
        cells = np.hstack((Qd,Td))
        self._cells = cells
        self._types = types
        self._F = F

    def _make_vertices(self):
        Fx = self._sides
        N = self.vectors.shape[0]
        Z = utilities.normalize(self.vectors)
        L = np.linalg.norm(self.vectors, axis=1)
        if self._position == 'tail':
            anchor = self.anchor+self.vectors*self._scale_factor+self._offset*Z
        elif self._position == 'head':
            anchor = self.anchor - self._offset*Z
        elif self._position == 'center':
            anchor = self.anchor+0.5*self.vectors*self._scale_factor+self._offset*Z
        if self._tip_length is None:
            tip_length = min([0.4*np.mean(L),0.6*np.min(L)])
        else:
            tip_length = self._tip_length
        if self._tip_radius is None:
            tip_radius = 0.15*tip_length
        else:
            tip_radius = self._tip_radius
        if self._shaft_radius is None:
            shaft_radius = 0.4*tip_radius
        else:
            shaft_radius = self._shaft_radius
        phi = np.linspace(0, 2*np.pi, Fx+1)[:-1]
        P0 = np.zeros((Fx,3))
        P1 = np.array([tip_radius*np.sin(phi), tip_radius*np.cos(phi),
                       -tip_length*np.ones(Fx)]).T
        P2 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                       -tip_length*np.ones(Fx)]).T
        P5 = np.array([0.3*tip_radius*np.sin(phi),0.3*tip_radius*np.cos(phi),
                       -0.3*tip_length*np.ones(Fx)]).T
        P6 = np.array([0.6*tip_radius*np.sin(phi),0.6*tip_radius*np.cos(phi),
                       -0.6*tip_length*np.ones(Fx)]).T
        P7 = np.array([0.1*tip_radius*np.sin(phi),0.1*tip_radius*np.cos(phi),
                       -0.1*tip_length*np.ones(Fx)]).T
        P3 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                       -np.ones(Fx)]).T
        P4 = np.array([[0,0,-1]])
        P = np.vstack((P0,P7,P5,P6,P1,P2,P3,P4)) * self._scale_factor
        Pt = np.zeros((N, P.shape[0],3))
        Pt[:] = P
        Off = np.tile(anchor, P.shape[0])
        Off = np.reshape(Off,(P.shape[0]*N, 3))
        Ld = np.repeat(L, Fx+1)
        Ld = np.reshape(Ld, (N, Fx+1))
        Pt[:,6*Fx:,2] *= Ld
        X = utilities.orthogonal_vectors(Z)
        Y = np.cross(Z, X, axis=1)
        J = np.zeros((N, 3, 3))
        J[:,0,:] = X
        J[:,1,:] = Y
        J[:,2,:] = Z
        Pt = np.einsum('ijk,ikl->ijl', Pt, J)
        Pt = np.reshape(Pt,(P.shape[0]*N, 3))
        points = Pt + Off
        self._lengths = L
        self._G = int(points.shape[0]//self.vectors.shape[0])
        self.vertices = points

    def cell_arrays(self):
        return self._cells, self._types


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                  Spheres
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class Spheres(object):

    def __init__(self, center, radius, **kwargs):

        self.center = center

        self.radius = radius

        self._resolution = kwargs.get('resolution', 20)

        self._scale_factor = kwargs.get('scale_factor', 1)

        self._mask = None

        self._make_faces()

        self._make_vertices()

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def F(self):
        return self._F

    @property
    def G(self):
        return self._G

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius):
        if type(radius) is int or type(radius) is float:
            radius = np.tile(radius, self.center.shape[0])
        try:
            radius.shape
        except AttributeError:
            radius = np.array(radius)
        self._radius = radius

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, center):
        center = np.array(center)
        self._center = center

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, resolution):
        if self._resolution != resolution:
            self._resolution = resolution
            self.make_faces()

    @property
    def mask(self):
        if self._mask is None:
            norm = np.linalg.norm(self._vectors, axis=1)
            self._mask = norm != 0
        return self._mask

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        update_faces = False
        if update_faces:
            self._make_faces()
        self._make_vertices()

    def _make_faces(self):
        Fx = 2*self._resolution; Fy = self._resolution-1
        N = self.center.shape[0]
        F = N * (Fx * (Fy + 1))
        M = np.arange(Fx*Fy, dtype=int)
        M = np.reshape(M,(Fy,Fx))
        M = np.insert(M,Fx,M[:,0],axis=1)
        Q = np.zeros((Fy-1, Fx, 4) ,dtype=int)
        Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
        Q[:,:,2] = M[:M.shape[0]-1,1:]
        Q[:,:,1] = M[1:,1:]
        Q[:,:,0] = M[1:,:M.shape[1]-1]
        T = np.zeros((2, Fx, 3) ,dtype=int)
        T[0,:,0] = M[0,:M.shape[1]-1] + 1
        T[0,:,1] = M[0,1:] + 1
        T[1,:,0] = M[-1,1:] + 1
        T[1,:,1] = M[-1,:M.shape[1]-1] + 1
        T[1,:,2] = Fx*Fy + 1
        O = np.arange(0, N*(Fx*Fy+2), Fx*Fy+2)
        O = np.reshape(O, (N,1,1))
        Q = Q.reshape(((Fx)*(Fy-1), 4),order='C') + 1
        T = T.reshape((2*Fx,3), order='C')
        t1 = np.repeat(1, N*Fx*(Fy-1))
        t2 = np.repeat(0, N*Fx)
        types = np.hstack((t1,t2,t2))
        Qd = np.zeros((N,Q.shape[0],Q.shape[1]),dtype=int)
        Qd[:] = Q
        Qd += O
        Td = np.zeros((N,T.shape[0],T.shape[1]),dtype=int)
        Td[:] = T
        Td += O
        Qd = np.insert(Qd,0,4,axis=2)
        Td = np.insert(Td,0,3,axis=2)
        Qd = np.reshape(Qd,(Qd.shape[0]*Qd.shape[1]*Qd.shape[2]))
        Td = np.reshape(Td,(Td.shape[0]*Td.shape[1]*Td.shape[2]))
        cells = np.hstack((Qd,Td))
        self._cells = cells
        self._types = types
        self._F = F

    def _make_vertices(self):
        Fx = 2*self._resolution; Fy = self._resolution-1
        N = self.center.shape[0]
        C = self.center
        phi = np.linspace(0, 2*np.pi, Fx+1)[:-1]
        theta = np.linspace(0, np.pi, Fy+2)
        theta = theta[1:theta.shape[0]-1]
        theta, phi, radius = np.meshgrid(theta, phi, self.radius)
        Px = radius*np.cos(phi)*np.sin(theta)
        Py = radius*np.sin(phi)*np.sin(theta)
        Pz = radius*np.cos(theta)
        Px = np.reshape(Px,(Fx*Fy*N), order='F') + np.repeat(C[:,0], Fx*Fy)
        Py = np.reshape(Py,(Fx*Fy*N), order='F') + np.repeat(C[:,1], Fx*Fy)
        Pz = np.reshape(Pz,(Fx*Fy*N), order='F') + np.repeat(C[:,2], Fx*Fy)
        i0 = np.arange(0, Px.shape[0], Fx*Fy)
        Px = np.insert(Px, i0, C[:,0])
        Py = np.insert(Py, i0, C[:,1])
        Pz = np.insert(Pz, i0, C[:,2] + self.radius)
        i0 = np.arange(Fx*Fy+1, Px.shape[0], Fx*Fy+1)
        i0 = np.hstack((i0, Px.shape[0]))
        Px = np.insert(Px, i0, C[:,0])
        Py = np.insert(Py, i0, C[:,1])
        Pz = np.insert(Pz, i0, C[:,2] - self.radius)
        P = np.vstack((Px,Py,Pz)).T
        self._G = int(P.shape[0]//self.center.shape[0])
        self.vertices = P

    def cell_arrays(self):
        return self._cells, self._types


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                          Mesh Polyline Pipe Glyph
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Pipe(object):

    def __init__(self, polyline, **kwargs):

        self._polyline = polyline

        self._radius = kwargs.get('radius', 0.1)

        self._sides = kwargs.get('sides', 20)

        self._comb = kwargs.get('comb', False)

        self._untwist = kwargs.get('untwist', False)

        self._make_faces()

        self._make_vertices()

    @property
    def type(self):
        return 'Glyph'

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def F(self):
        return self._types.shape[0]

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius):
        if self._radius != radius:
            self._radius = radius
            self._make_vertices()

    @property
    def sides(self):
        return self._sides

    @sides.setter
    def sides(self, sides):
        if self._sides != sides:
            self._sides = sides
            self._make_faces()
            self._make_vertices()

    @property
    def polyline(self):
        return self._polyline

    @polyline.setter
    def polyline(self, polyline):
        self._polyline = polyline
        self._make_faces()
        self._make_vertices()

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        update_faces = False
        polyline = kwargs.get('polyline', self._polyline)
        if self._polyline.V != polyline.V:
            update_faces = True
        if self._polyline.closed != polyline.closed:
            update_faces = True
        self._polyline = polyline
        self._radius = kwargs.get('radius', self._radius)
        if self._sides != kwargs.get('sides', self._sides):
            self._sides = kwargs.get('sides', self._sides)
            update_faces = True
        if update_faces:
            self._make_faces()
        self._make_vertices()

    def _make_faces(self):
        Fx = self._sides
        Fy = self._polyline.V
        if not self._polyline.closed:
            Fy = self._polyline.V
            M = np.arange(Fx*Fy, dtype=int)
            M = np.reshape(M,(Fy,Fx))
            M = np.insert(M,Fx,M[:,0],axis=1)
            Q = np.zeros((Fy-1, Fx, 4) ,dtype=int)
            Q[:,:,0] = M[:M.shape[0]-1,:M.shape[1]-1]
            Q[:,:,1] = M[:M.shape[0]-1,1:]
            Q[:,:,2] = M[1:,1:]
            Q[:,:,3] = M[1:,:M.shape[1]-1]
            Q = Q.reshape(((Fx)*(Fy-1),4),order='C')
        else:
            M = np.arange(Fx*Fy, dtype=int)
            M = np.reshape(M,(Fy,Fx))
            M = np.insert(M,Fx,M[:,0],axis=1)
            M = np.insert(M,Fy,M[0,:],axis=0)
            Q = np.zeros((Fy, Fx, 4) ,dtype=int)
            Q[:,:,0] = M[:M.shape[0]-1,:M.shape[1]-1]
            Q[:,:,1] = M[:M.shape[0]-1,1:]
            Q[:,:,2] = M[1:,1:]
            Q[:,:,3] = M[1:,:M.shape[1]-1]
            Q = Q.reshape(((Fx)*(Fy),4),order='C')
        types = np.ones(Q.shape[0], dtype=int)
        lengths = 4 * np.ones((Q.shape[0],1))
        cells = np.hstack((lengths, Q))
        self._cells = np.array(cells.flatten('C'), 'i')
        self._types = types

    def _make_vertices(self):
        Fx = self._sides
        Fy = self._polyline.V
        V1, C = self._polyline.vertex_bisectors(return_cosines=True)
        if self._comb:
            V1 = comb_normals(self._polyline)
        if self._untwist:
            Va = V1[1:self.polyline.V]
            Vb = V1[0:self.polyline.V-1]
            D = np.einsum('ij,ij->i', Va, Vb)
            sign = 1
            for i in range(len(D)):
                V1[i] *= sign
                if D[i] < -0.2:
                    sign *= -1
            V1[-1] = V1[-2]
        V2 = self._polyline.vertex_tangents()
        V3 = np.cross(V1, V2)
        B = np.sqrt(np.abs(0.5*(1-C))) + 1e-10
        B = self._radius / B
        B = np.hstack((B,B,B))
        phi = np.linspace(0,2*np.pi,Fx+1)
        phi = np.delete(phi,-1)
        r = self._radius
        i = np.repeat(np.arange(Fy), Fx)
        phi = np.tile(phi, Fy)
        Vi = self._polyline.vertices
        Px = B[i]*np.cos(phi)*V1[i,0] + r*np.sin(phi)*V3[i,0] + Vi[i,0]
        Py = B[i]*np.cos(phi)*V1[i,1] + r*np.sin(phi)*V3[i,1] + Vi[i,1]
        Pz = B[i]*np.cos(phi)*V1[i,2] + r*np.sin(phi)*V3[i,2] + Vi[i,2]
        Px = np.reshape(Px,(Fx*(Fy)), order='C')
        Py = np.reshape(Py,(Fx*(Fy)), order='C')
        Pz = np.reshape(Pz,(Fx*(Fy)), order='C')
        points = np.vstack((Px,Py,Pz)).T
        self.vertices = points

    def cell_arrays(self):
        return self._cells, self._types

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                          Mesh Circle Disc
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Discs(object):

    def __init__(self, circle, **kwargs):

        self._circle = circle

        self._cells = None

        self._types = None

        self._make_vertices()

        self._make_faces()



    @property
    def type(self):
        return 'Glyph'

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def F(self):
        return self._types.shape[0]

    @property
    def vertices(self):
        return self._circle.vertices

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        pass

    def _make_faces(self):
        F = np.arange(self.V)
        N = len(self._circle.frame.origin)
        F = np.reshape(F, (N, self._circle.sampling))
        lengths = self._circle.sampling * np.ones(N)
        cells = np.column_stack((lengths, F))
        self._cells = np.array(cells.flatten('C'), 'i')
        types = 2*np.ones(N, dtype=int)
        self._types = types

    def _make_vertices(self):
        self._circle.make_vertices()

    def cell_arrays(self):
        return self._cells, self._types


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                          Mesh Circle Disc
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Cones(object):

    def __init__(self, circle, height=[1], truncate=0):

        self._circle = circle

        self._height = height

        self._cells = None

        self._types = None

        self._truncate = truncate

        self._make_faces()

        self._make_vertices()

    @property
    def type(self):
        return 'Glyph'

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def F(self):
        return self._types.shape[0]

    @property
    def vertices(self):
        return self._vertices

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        pass

    def _make_faces(self):
        N = self._circle.frame.origin.shape[0]
        S = self._circle.sampling
        F = np.arange(N*S)
        F = F.reshape((N,S))
        T = N*S + F
        F = np.dstack((F, np.roll(F,1, axis=1), np.roll(T,1, axis=1), T))
        F = F.reshape(((N*S,4)))
        F = np.column_stack((np.repeat(4,N*S), F))
        self._cells = np.array(F.flatten('C'), 'i')
        self._types = np.repeat(1, N*S)

    def _make_vertices(self):
        self._circle.make_vertices()
        H = np.column_stack((self._height,self._height,self._height))
        A = self._circle.frame.origin + H*self._circle.frame.e3
        A = np.repeat(A, self._circle.sampling, axis=0)
        A = (1-self._truncate)*A + (self._truncate)*self._circle.vertices
        self._vertices = np.vstack((self._circle.vertices, A))

    def cell_arrays(self):
        return self._cells, self._types


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                          Mesh Circle Disc
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Rings(object):

    def __init__(self, circle, **kwargs):

        self._circle = circle

        self._radius = kwargs.get('radius', 0.1)

        self._sides = kwargs.get('sides', 20)

        self._cells = None

        self._types = None

        self._vertices = None

        self._make_faces()

        self._make_vertices()

    @property
    def type(self):
        return 'Glyph'

    @property
    def V(self):
        return self.vertices.shape[0]

    @property
    def F(self):
        return self._types.shape[0]

    @property
    def vertices(self):
        return self._vertices

    def update(self, **kwargs):
        pass

    def _make_faces(self):
        Fy = self._circle.sampling
        Fx = self._sides
        N = self._circle.frame.origin.shape[0]
        F = np.arange(N*Fy*Fx)
        F = F.reshape((N,Fx,Fy))
        Frx = np.roll(F,-1,axis=1)
        F = F.reshape((Fx*N,Fy))
        Frx = Frx.reshape((Fx*N,Fy))
        Q = np.dstack((F, np.roll(F,1,axis=1), np.roll(Frx,1,axis=1), Frx))
        Q = np.column_stack((4*np.ones(N*Fx*Fy, 'i'), Q.reshape((N*Fx*Fy,4))))
        self._cells = Q.flatten('C')
        self._types = np.ones(N*Fx*Fy, 'i')
        self._F = N*Fx*Fy

    def _make_vertices(self):
        r = self._radius
        Fy = self._circle.sampling
        Fx = self._sides
        N = self._circle.frame.origin.shape[0]
        C = self._circle.frame.origin
        u = np.linspace(0, 2*np.pi, Fx+1)[:-1]
        v = np.linspace(0, 2*np.pi, Fy+1)[:-1]
        u, R, v = np.meshgrid(u, self._circle.radius, v)
        r = r*np.ones(R.shape)
        Px = r*np.cos(u)*np.cos(v) + R*np.cos(v)
        Py = r*np.cos(u)*np.sin(v) + R*np.sin(v)
        Pz = r*np.sin(u)
        Px = np.reshape(Px,(Fx*Fy*N), order='C')
        Py = np.reshape(Py,(Fx*Fy*N), order='C')
        Pz = np.reshape(Pz,(Fx*Fy*N), order='C')
        X = np.repeat(self._circle.frame.e1, Fx*Fy, axis=0)
        Y = np.repeat(self._circle.frame.e2, Fx*Fy, axis=0)
        Z = np.repeat(self._circle.frame.e3, Fx*Fy, axis=0)
        Px = np.column_stack((Px,Px,Px))
        Py = np.column_stack((Py,Py,Py))
        Pz = np.column_stack((Pz,Pz,Pz))
        P = Px*X + Py*Y + Pz*Z + np.repeat(C, Fx*Fy, axis=0)
        self._G = int(P.shape[0]//N)
        self._vertices = P

    def cell_arrays(self):
        return self._cells, self._types


