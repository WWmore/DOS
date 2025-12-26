#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.geometry.frame import Frame

# -----------------------------------------------------------------------------

'''circle.py: The Circle data structure class'''

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Circle
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class Circle(object):

    def __init__(self, frame=Frame(), radius=[1], sampling=50):

        self.name = 'circle'

        self.frame = frame

        self.radius = radius

        self.sampling = int(sampling)

        self.vertices = None

        self._type = 'Circle'

    @property
    def type(self):
        return self._type

    @property
    def V(self):
        if self.vertices is not None:
            return len(self.vertices)
        else:
            return None

    @property
    def E(self):
        if self.vertices is not None:
            return len(self.vertices)
        else:
            return None

    @property
    def radius(self):
        try:
            return self._radius
        except AttributeError:
            return None

    @radius.setter
    def radius(self, radius):
        if type(radius) is int or type(radius) is float:
            radius = [radius]
        self._radius = np.array(radius, dtype=float)

    # -------------------------------------------------------------------------
    #                               Sampling
    # -------------------------------------------------------------------------

    def make_vertices(self):
        N = self.frame.origin.shape[0]
        phi = np.linspace(0, 2*np.pi - 2*np.pi / self.sampling, self.sampling)
        phi = np.tile(phi, N)
        r = np.repeat(self.radius, self.sampling)
        Ox = np.repeat(self.frame.origin[:, 0], self.sampling)
        Oy = np.repeat(self.frame.origin[:, 1], self.sampling)
        Oz = np.repeat(self.frame.origin[:, 2], self.sampling)
        e1x = np.repeat(self.frame.e1[:, 0], self.sampling)
        e1y = np.repeat(self.frame.e1[:, 1], self.sampling)
        e1z = np.repeat(self.frame.e1[:, 2], self.sampling)
        e2x = np.repeat(self.frame.e2[:, 0], self.sampling)
        e2y = np.repeat(self.frame.e2[:, 1], self.sampling)
        e2z = np.repeat(self.frame.e2[:, 2], self.sampling)
        vx = Ox + r*np.sin(phi)*e1x + r*np.cos(phi)*e2x
        vy = Oy + r*np.sin(phi)*e1y + r*np.cos(phi)*e2y
        vz = Oz + r*np.sin(phi)*e1z + r*np.cos(phi)*e2z
        vertices = np.array([vx, vy, vz]).T
        self.vertices = vertices
        self._type = 'Polyline'

    def cell_array(self):
        vi = np.arange(self.V)
        vj = np.roll(vi, -1)
        c = np.repeat(2, self.vertices.shape[0])
        i = np.arange(0, self.V, self.sampling)
        j = np.arange(self.sampling, self.V + 1, self.sampling)
        j = j - 1
        vj[j] = vi[i]
        cells = np.vstack((c, vi, vj)).T
        cells = np.ravel(cells)
        return cells


# -----------------------------------------------------------------------------
#                                 Generation
# -----------------------------------------------------------------------------


def circle_three_points(p1, p2, p3, center=False): #Hui: add center
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    t = p2-p1
    u = p3-p1
    v = p3-p2
    n = np.cross(t, u)
    nsl = np.linalg.norm(n, axis=1)**2 #Huinote: problem for only three points
    nsl[np.where(nsl == 0)] == 1e-20
    insl2 = 1.0 / (2.0*nsl)
    insl3 = np.array([insl2, insl2, insl2]).T
    a = np.einsum('ij,ij->i', t, t) * np.einsum('ij,ij->i', u, v)
    b = np.einsum('ij,ij->i', t, v) * np.einsum('ij,ij->i', u, u)
    c = np.array([a, a, a]).T * u - np.array([b, b, b]).T * t
    C = p1 + insl3 * c
    r = (np.einsum('ij,ij->i', t, t) * np.einsum('ij,ij->i', u, u) *
         np.einsum('ij,ij->i', v, v) * insl2*0.5)**0.5
    n = n/np.linalg.norm(n, axis=1, keepdims=True)
    e2 = t/np.linalg.norm(t, axis=1, keepdims=True)
    e1 = np.cross(e2, n)
    frame = Frame(C, e1, e2, n)
    if center:
        return Circle(frame,r), C, r
    else:
        return Circle(frame, r)