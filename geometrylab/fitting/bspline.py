# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# Create a 3-dimensional B-spline Curve


# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'

class BSpline(object):

    def __init__(self, degree=3, control_points=None, closed=False):

        self.degree = degree

        self.control_points = control_points

        self.closed = closed

        self.knot_vector = None

        if control_points is not None:
            self.make_knot_vector()

    def basis_functions(self, t, **kwargs):
        control_points = kwargs.get('control_points', self.control_points)
        n = kwargs.get('degree', self.degree)
        kv = kwargs.get('knot_vector', self.knot_vector)
        m = len(control_points) - 1
        N = np.zeros((len(t), m + 2))
        for i in range(m + 1):
            c1 = t >= kv[i]
            c2 = t < kv[i+1]
            ind = np.logical_and(c1, c2)
            N[ind,i] = 1
        for r in range(1, n + 1):
            for j in range(m + 1):
                A = (t - kv[j]) / (kv[j+r] - kv[j] + 1e-10)
                B = (kv[j+r+1] - t) / (kv[j+r+1] - kv[j+1] + 1e-10)
                N[:,j] = A*N[:,j] + B*N[:,j+1]
        N = np.delete(N, -1, axis=1)
        return N

    def derivatives(self, t):
        P = np.copy(self.control_points)
        knot_vector = self.knot_vector
        n = self.degree
        D = []
        for i in range(2):
            P0 = np.delete(P, -1, axis=0)
            P1 = np.delete(P, 0, axis=0)
            p = len(P) - 1
            kv = np.array(knot_vector)
            P = np.einsum('i,ij->ij', n/(kv[1+n:1+p+n] - kv[1:p+1]), P1 - P0)
            n = n - 1
            knot_vector = knot_vector[1:-1]
            N = self.basis_functions(t, degree=n,
                                            knot_vector=knot_vector,
                                            control_points=P)
            D.append(np.dot(N, P))
        return D[0], D[1]

    def frame(self, t):
        F, S = self.derivatives(t)
        A = np.einsum('ij,ij->i', F, S) / np.linalg.norm(F, axis=1)**2
        N = S - np.einsum('i,ij->ij', A, F)
        N = -N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-10)
        T = F / np.linalg.norm(F, axis=1, keepdims=True)
        return T, N

    def points(self, t):
        N = self.basis_functions(t)
        points = np.dot(N, self.control_points)
        return points

    def make_knot_vector(self, t_min=0, t_max=1):
        m = len(self.control_points) - 1
        n = self.degree
        kv = []
        for i in range(n + 1):
            kv.append(t_min)
        for i in range(m - n):
            step = float(t_max - t_min) / (m - n + 1)
            kv.append(t_min + step*(i + 1))
        for i in range(n + 1):
            kv.append(t_max)
        self.knot_vector = kv



import os, sys
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#print(path)
sys.path.append(path)

from geometrylab import vtkplot




P = np.array([[0,0,0], [0,10,0], [10,-10,0], [10,10,10], [-10,-10,0]])
sp = BSpline(control_points=P, degree=3)
t = np.array([0.0001, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99999])
t = np.random.random(200)

Pt = sp.points(t)

T, N = sp.frame(t)
pl_v = vtkplot.Vectors(T,anchor=Pt)

pl_p = vtkplot.Points(P, radius=0.3, color='r')
pl_pt = vtkplot.Points(Pt, radius=0.01)


vtkplot.view([pl_p,pl_pt, pl_v])













