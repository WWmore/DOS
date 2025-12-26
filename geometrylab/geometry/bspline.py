# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import spatial

#------------------------------------------------------------------------------

from geometrylab.fitting.linearregression import linear_regression

from geometrylab.geometry import Polyline

#------------------------------------------------------------------------------




__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------


class BSpline(object):

    def __init__(self, degree=3, control_points=None, closed=False,
                 knot_vector=None):

        self.name = 'bspline'

        self._degree = degree

        self._kdtree = None

        self._control_points = control_points

        self._closed = closed

        self._knot_vector = knot_vector

        self._epsilon = 5e-4

        self.sampling = 200

        if control_points is not None and knot_vector is None:
            self.make_knot_vector()

    def __str__(self):
        name = 'B-spline curve: '
        info = 'degree = {}, number of points = {}, dimension = {}'
        out = name + info.format(self.n, self.m+1, self.d)
        return out


    @property
    def type(self):
        return 'BSpline'

    @property
    def vertices(self):
        return self.control_points

    @property
    def degree(self):
        return self._degree

    @degree.setter
    def degree(self, degree):
        self._degree = degree
        self._kdtree = None

    @property
    def closed(self):
        return self._closed

    @closed.setter
    def closed(self, bool):
        self._closed = bool
        self._kdtree = None

    @property
    def knot_vector(self):
        return self._knot_vector

    @knot_vector.setter
    def knot_vector(self, knot_vector):
        self._knot_vector = knot_vector
        self._kdtree = None

    @property
    def control_points(self):
        return self._control_points

    @control_points.setter
    def control_points(self, control_points):
        self._control_points = control_points
        self._kdtree = None

    @property
    def epsilon(self):
        return self._epsilon

    @epsilon.setter
    def epsilon(self, epsilon):
        self._epsilon = epsilon
        self._kdtree = None

    @property
    def m(self):
        return len(self.control_points) - 1

    @property
    def n(self):
        return int(self.degree)

    @property
    def d(self):
        return self.control_points.shape[1]

    @property
    def t_min(self):
        return min(self.knot_vector)

    @property
    def t_max(self):
        return max(self.knot_vector)

    def basis_functions(self, t, **kwargs):
        n = kwargs.get('degree', self.degree)
        kv = kwargs.get('knot_vector', self.knot_vector)
        m = kwargs.get('m', self.m)
        N = np.zeros((len(t), m + 2))
        for i in range(m + 1):
            c1 = t >= kv[i]
            c2 = t < kv[i+1]
            ind = np.logical_and(c1, c2)
            N[ind,i] = 1
            N[t==kv[-1],m] = 1
        for r in range(1, n + 1):
            for j in range(m + 1):
                a = kv[j+r] - kv[j]
                if a != 0:
                    A = (t - kv[j]) / a
                else:
                    A = np.zeros(len(t))
                b = kv[j+r+1] - kv[j+1]
                if b != 0:
                    B = (kv[j+r+1] - t) / b
                else:
                    B = np.zeros(len(t))
                N[:,j] = A*N[:,j] + B*N[:,j+1]
        N = np.delete(N, -1, axis=1)
        return N

    def first_derivative_basis_functions(self, t):
        d = 1
        kv = np.array(self.knot_vector)
        N = self.basis_functions(t, degree = self.n - d,
                                knot_vector = self.knot_vector[d:-d],
                                m = self.m-d)
        i = np.arange(self.m+1)
        a = (kv[i+self.n] - kv[i])
        zero = np.where(a == 0)[0]
        a[zero] = 1
        f = (self.n-d+1) / a
        f[zero] = 0
        f = np.hstack((f,[0]))
        N = np.column_stack((np.zeros(len(N)), N, np.zeros(len(N))))
        N = np.einsum('j,ij->ij', f, N)
        D = N[:,0:-1] - N[:,1:]
        return D

    def stretch_energy_matrix(self, resolution=1000):
        t = np.linspace(self.t_min, self.t_max, resolution)
        B = self.first_derivative_basis_functions(t)
        B = np.einsum('ki,kj->ijk', B, B)
        B = np.sum(B, axis=2)
        return B

    def hodograph(self):
        P = np.copy(self.control_points)
        n = self.degree
        P0 = np.delete(P, -1, axis=0)
        P1 = np.delete(P, 0, axis=0)
        kv = np.array(self.knot_vector)
        m = self.m
        P = np.einsum('i,ij->ij', n/(kv[1+n:1+m+n] - kv[1:m+1]), P1 - P0)
        H = BSpline(degree=n-1, control_points=P,
                    knot_vector=self.knot_vector[1:-1])
        return H

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
            N = self.basis_functions(t, degree = n,
                                    knot_vector = knot_vector,
                                    m = len(P)-1)
            D.append(np.dot(N, P))
        return D[0], D[1]

    def frame(self, t):
        F, S = self.derivatives(t)
        T =  F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-8)
        if self.d == 2:
            N = np.column_stack((T[:,1], -T[:,0]))
        elif self.d == 3:
            A = np.einsum('ij,ij->i', F, S) / np.linalg.norm(F, axis=1)**2
            N = S - np.einsum('i,ij->ij', A, F)
            N = -N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-8)
        return T, N

    def points(self, t):
        N = self.basis_functions(t)
        points = np.dot(N, self.control_points)
        return points

    def make_knot_vector(self, t_min=0, t_max=1):
        n = self.degree
        kv = []
        for i in range(n + 1):
            kv.append(t_min)
        for i in range(self.m - n):
            step = float(t_max - t_min) / (self.m - n + 1)
            kv.append(t_min + step*(i + 1))
        for i in range(n + 1):
            kv.append(t_max)
        self.knot_vector = kv

    def point_curve_distance(self, points):
        t = self.closest_points(points)
        distance = np.linalg.norm(self.points(t) - points, axis=1)
        return distance

    def closest_points(self, points):
        n = int(1./self._epsilon)
        t = np.linspace(self.t_min, self.t_max, n)
        curve_points = self.points(t)
        if self._kdtree is None:
            self._kdtree = spatial.cKDTree(curve_points)
        closest = self._kdtree.query(points, k=2)[1]
        i0 = np.minimum(closest[:,0], closest[:,1])
        i1 = np.maximum(closest[:,0], closest[:,1])
        P0 = curve_points[i0]
        P1 = curve_points[i1]
        V = (P1 - P0) / np.linalg.norm(P1 - P0, axis=1,  keepdims=True)
        dt = np.einsum('ij,ij->i', points - P0, V)
        tc = t[i0] + dt * ((self.t_max - self.t_min) / n)
        tc[tc > self.t_max] = self.t_max
        tc[tc < self.t_min] = self.t_min
        return tc

    def sample_points(self, n=None):
        if n is None:
            n = self.sampling
        t = np.linspace(min(self.knot_vector), max(self.knot_vector), n)
        return self.points(t)

    def fit_to_points(self, points, iterations=10, straightening=0.1):
        for i in range(iterations):
            t = self.closest_points(points)
            B = self.basis_functions(t)
            T, N = self.frame(t)
            P = np.zeros((self.m+1, 0))
            for j in range(self.d):
                M = np.einsum('i,ij->ij', N[:,j], B)
                R = -2 * np.eye(self.m+1)
                R[0,0] = 0
                R[-1,-1] = 0
                k = np.arange(1,self.m)
                R[k,k-1] = 1
                R[k,k+1] = 1
                a = points[:,j] * N[:,j]
                M = np.vstack((M, straightening*R))
                a = np.hstack((a, np.zeros(self.m+1)))
                Pi = np.linalg.lstsq(M, a, rcond=-1)[0]
                P = np.column_stack((P,Pi))
            self.control_points = P


    def fit_line_to_points(self, points, control_points_number=None, delta=0.1):
        if control_points_number is None:
            control_points_number = self.n + 1
        N = control_points_number
        regression = linear_regression(points)
        L = regression[0][-1]
        C = regression[1]
        t = np.einsum('j,ij->i', L, points - C)
        e = (np.max(t) - np.min(t)) * delta
        u = np.linspace(np.min(t) - e, np.max(t) + e, N)
        P = np.einsum('i,j->ij', u, L)
        self.control_points = P[:,:] + C[:]
        self.make_knot_vector()

    def curve_polyline(self):
        return Polyline(self.sample_points(), closed=self.closed)

    def control_polyline(self):
        return Polyline(self.control_points, closed=self.closed)





#------------------------------------------------------------------------------
if __name__ == '__main__':

    from geometrylab import vtkplot
    
    if 1:
        "plot Bezier curve, control points and control polygon:"
        P = np.array([[3,0,0], [0,10,0], [10,-10,0], [10,10,0]])
        sp = BSpline(control_points=P, degree=3)
        #sp = sp.hodograph()
        #sp = sp.hodograph()
        
        crv = sp.curve_polyline()
        ctrl = sp.control_polyline()
        
        pl_crv = vtkplot.Polyline(polyline=crv, tube_radius=0.1, color='r',sampling=500)
        pl_ctrl = vtkplot.Polyline(polyline=ctrl, color='black')
        pl_pts = vtkplot.Points(P, radius=0.5, color='black')
    
        vtkplot.view([pl_pts, pl_crv, pl_ctrl])
    
    else:
        "plot B-spline curve, control points and control polygon: "
        P = np.array([[3,0,0], [0,10,0], [10,-10,0], [10,10,0], [-10,-10,0]])
        sp = BSpline(control_points=P, degree=3,closed=False)

        crv = sp.curve_polyline()
        ctrl = sp.control_polyline()
        
        pl_crv = vtkplot.Polyline(polyline=crv, tube_radius=0.1, color='r',sampling=500)
        pl_ctrl = vtkplot.Polyline(polyline=ctrl, color='black')
        pl_pts = vtkplot.Points(P, radius=0.5, color='black')
    
        vtkplot.view([pl_pts, pl_crv, pl_ctrl])
