#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import sparse

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.geometry import Frame

from geometrylab.geometry import Circle

from geometrylab.optimization.guidedprojectionbase import GuidedProjectionBase

# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'

'''
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
                          SELF SUPPORTING OPTIMIZATION
-------------------------------------------------------------------------------
                'Form-finding with Polyhedral Meshes Made Simple'
                             (Tang et al., 2014)
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
                                     Notes
-------------------------------------------------------------------------------
 indexing:
 V = vertices
 E = edges
 F = faces
-------------------------------------------------------------------------------
 unknows vector:
 3V | 3F | 3E | 4F | 3F

                | 3V              | N1             | N2                  | N3
 [x_V, y_V, z_V, nx_F, ny_F, nz_F, l_E, w_E, lam_E, Ax_F, Ay_F, Az_F, A_F,
 Cx_F, Cy_F, Cz_F]
-------------------------------------------------------------------------------
 sign convenction: w_ij > 0 = compression
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
                                  Constraints
-------------------------------------------------------------------------------
                      fi = Hi_jk X_j X_k + bi_j X_j - ci
-------------------------------------------------------------------------------
            H_ij = Hi_jk X_k + bi_j;    r_i = 1/2 Hi_jk X_j X_k - ci
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
'''

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class GuidedProjection(GuidedProjectionBase):

    _N1 = 0

    _N2 = 0

    _N3 = 0

    _N4 = 0

    _compression = 0

    _mesh = None

    def __init__(self):
        GuidedProjectionBase.__init__(self)

        weights = {

        'planarity' : 0,

        'normal' : 0,

        'edge_length' : 0,

        'area' : 0,

        'geometric' : 1,

        'circularity' : 0,

        'fixed_vertices' : 1,

        'fixed_corners' : 0,

        'gliding' : 0, # Huinote: glide on boundary, used for itself boundary

        'reference_closeness' : 0,

        'self_closeness' : 0,

        'boundary_closeness' : 0,

        'equilibrium' : 0,

        'fixed_boundary_normals': 0,

        'mesh_fairness' : 0.01,

        'tangential_fairness' : 0.01,

        'boundary_fairness' : 0.01,

        'spring_fairness' : 0}

        self.add_weights(weights)

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh
        self.initialization()

    @property
    def compression(self):
        return self._compression

    @compression.setter
    def compression(self, bool):
        if bool:
            self._compression = 2*self.get_weight('equilibrium')
            self.reinitialize = True
        elif self._compression > 0:
             self._compression = 0

    @property
    def tension(self):
        return -self._compression

    @tension.setter
    def tension(self, bool):
        if bool:
            self._compression = -2*self.get_weight('equilibrium')
            self.reinitialize = True
        elif self._compression < 0:
            self._compression = 0

    @property
    def max_weight(self):
        return max(self.get_weight('equilibrium'),
                   self.get_weight('planarity'),
                   self.get_weight('geometric'))

    #
    #--------------------------------------------------------------------------
    #                               Initialization
    #--------------------------------------------------------------------------

    def set_weights(self):
        if self.get_weight('equilibrium') != 0:
            if self.mesh.area_load != 0:
                if self.get_weight('area') == 0:
                    self.set_weight('area', 1 * self.get_weight('equilibrium'))
            else:
                self.set_weight('area', 0)
        if self.reinitialize:
            if self.get_weight('equilibrium') != 0:
                self.set_weight('edge_length', 1 * self.get_weight('equilibrium'))
                if self.mesh.area_load != 0:
                    self.set_weight('area', 1*self.get_weight('equilibrium'))
            if self.get_weight('planarity') != 0:
                self.set_weight('normal', 1*self.get_weight('planarity'))
        self.set_weight('fixed_vertices', 10 * self.max_weight)
        self.set_weight('gliding', 10 * self.max_weight)
        if self.get_weight('fixed_corners') != 0:
            self.set_weight('fixed_corners', 10 * self.max_weight)
        if self.get_weight('fixed_boundary_normals') != 0:
            self.set_weight('fixed_boundary_normals', 10 * self.max_weight)

    def set_dimensions(self):
        V = self.mesh.V
        F = self.mesh.F
        E = self.mesh.E
        N = 3*V
        N1 = 3*V
        N2 = 3*V
        N3 = 3*V
        N4 = 3*V
        if self.get_weight('planarity') != 0:
            N += 3*F
            N1 += 3*F
            N2 += 3*F
            N3 += 3*F
            N4 += 3*F
        if self.get_weight('equilibrium') != 0:
            N += 3*E
            N2 += 3*E
            N3 += 3*E
            N4 += 3*E
        if self.get_weight('area') != 0:
            N += 4*F
            N3 += 4*F
            N4 += 4*F
        if self.get_weight('circularity') != 0:
            N += 3*F
            N4 += 3*F
        if N1 != self._N1 or N2 != self._N2:
            self.reinitialize = True
        if N3 != self._N3 or N4 != self._N4:
            self.reinitialize = True
        if self._N2 - self._N1 == 0 and N2 - N1 > 0:
            self.mesh.reinitialize_densities()
        self._N = N
        self._N1 = N1
        self._N2 = N2
        self._N3 = N3
        self._N4 = N4

    def initialize_unknowns_vector(self):
        X = self.mesh.vertices.flatten('F')
        if self.get_weight('planarity') != 0:
            normals = self.mesh.face_normals()
            normals = normals.flatten('F')
            X = np.hstack((X, normals))
        if self.get_weight('equilibrium') != 0:
            lengths = self.mesh.edge_lengths()
            W = self.mesh.force_densities
            X = np.hstack((X, lengths, W, np.abs(W)**0.5))
        if self.get_weight('area') != 0:
            vector_area = self.mesh.face_vector_areas()
            face_area = np.linalg.norm(vector_area, axis=1)
            vector_area = vector_area.flatten('F')
            X = np.hstack((X, vector_area, face_area))
        if self.get_weight('circularity') != 0:
            centers = self.mesh.face_barycenters()
            centers = centers.flatten('F')
            X = np.hstack((X, centers))
            v, ej = self.mesh.vertex_ring_edges_iterators(sort=True)
        self._X = X
        self._X0 = np.copy(X)

    def make_errors(self):
        self.edge_length_error()
        self.equilibrium_error()
        self.face_areas_error()
        self.face_normals_error()
        self.face_vector_areas_error()
        self.planarity_error()
        self.geometric_error()

    def post_iteration_update(self):
        V = self.mesh.V
        E = self.mesh.E
        N1 = self._N1
        self.mesh.vertices[:,0] = self.X[0:V]
        self.mesh.vertices[:,1] = self.X[V:2*V]
        self.mesh.vertices[:,2] = self.X[2*V:3*V]
        if self.get_weight('equilibrium')!= 0:
            self.mesh.force_densities = self.X[N1+E:N1+2*E]
        else:
            self.mesh.force_densities = np.zeros(self.mesh.E)

    def on_reinitilize(self):
        self.mesh.reinitialize_force_densities()

    #--------------------------------------------------------------------------
    #                                  Results
    #--------------------------------------------------------------------------

    def vertices(self):
        V = self.mesh.V
        vertices = self.X[0:3*V]
        vertices = np.reshape(vertices, (V,3), order='F')
        return vertices

    def edge_lengths(self, initialized=False):
        if self.get_weight('equilibrium') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        E = self.mesh.E
        N1 = self._N1
        return X[N1:N1+E]

    def face_normals(self, initialized=False):
        if self.get_weight('planarity') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        V = self.mesh.V
        F = self.mesh.F
        normals = X[3*V:3*V+3*F]
        normals = np.reshape(normals, (F,3), order='F')
        return normals

    def face_vector_areas(self, initialized=False):
        if self.get_weight('area') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        F = self.mesh.F
        N2 = self._N2
        areas = X[N2:N2+3*F]
        areas = np.reshape(areas, (F,3), order='F')
        return areas

    def face_areas(self, initialized=False):
        if self.get_weight('area') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        F = self.mesh.F
        N2 = self._N2
        areas = X[N2+3*F:N2+4*F]
        return areas

    def force_densities(self):
        if self.get_weight('equilibrium') == 0:
            return None
        E = self.mesh.E
        N1 = self._N1
        return self.X[N1+E:N1+2*E]

    def face_circumcenters(self):
        if self.get_weight('circularity') == 0:
            return None
        X = self.X
        F = self.mesh.F
        O = self._N3
        centers = X[O:O+3*F]
        centers = np.reshape(centers, (F,3), order='F')
        return centers

    #--------------------------------------------------------------------------
    #                                  Errors
    #--------------------------------------------------------------------------

    def planarity_error(self):
        if self.get_weight('planarity') == 0:
            return None
        P = self.mesh.face_planarity()
        Emean = np.mean(P)
        Emax = np.max(P)
        self.add_error('planarity', Emean, Emax, self.get_weight('planarity'))

    def face_normals_error(self):
        if self.get_weight('planarity') == 0:
            return None
        N = self.face_normals()
        N0 = self.face_normals(initialized=True)
        norm = np.mean(np.linalg.norm(N, axis=1))
        Err = (np.linalg.norm(N-N0, axis=1)) / norm
        Emean = np.mean(Err)
        Emax = np.max(Err)
        self.add_error('face_normal', Emean, Emax, self.get_weight('normal'))
        return Emean, Emax

    def edge_length_error(self):
        if self.get_weight('equilibrium') == 0:
            return None
        L = self.edge_lengths()
        L0 = self.edge_lengths(initialized=True)
        norm = np.mean(L)
        Err = np.abs(L-L0) / norm
        Emean = np.mean(Err)
        Emax = np.max(Err)
        self.add_error('edge_length', Emean, Emax, self.get_weight('edge_length'))

    def face_vector_areas_error(self):
        if self.get_weight('area') == 0:
            return None
        A = self.face_vector_areas()
        A0 = self.face_vector_areas(initialized=True)
        norm = np.mean(np.linalg.norm(A, axis=1))
        Err = (np.linalg.norm(A-A0, axis=1)) / norm
        Emean = np.mean(Err)
        Emax = np.max(Err)
        self.add_error('face_vector_area', Emean, Emax, self.get_weight('area'))

    def face_areas_error(self):
        if self.get_weight('area') == 0:
            return None
        A = self.face_areas()
        A0 = self.face_areas(initialized=True)
        norm = np.mean(A)
        Err = (np.abs(A-A0)) / norm
        Emean = np.mean(Err)
        Emax = np.max(Err)
        self.add_error('face_area', Emean, Emax, self.get_weight('area'))

    def equilibrium_error(self):
        if self.get_weight('equilibrium') == 0:
            return None
        Err = self.mesh.equilibrium_error()
        Emean = np.mean(Err)
        Emax = np.max(Err)
        self.add_error('equilibrium', Emean, Emax, self.get_weight('equilibrium'))

    def geometric_error(self):
        if len(self._errors) == 0:
            return None
        n = 0
        geo_mean = 0
        geo_max = 0
        if self.get_weight('planarity') != 0:
            err = self.get_error('face_normal')
            geo_mean += err[0]
            geo_max = max([geo_max, err[1]])
            n += 1
        if self.get_weight('equilibrium') != 0:
            err = self.get_error('edge_length')
            geo_mean += err[0]
            geo_max = max([geo_max, err[1]])
            n += 1
            if self.get_weight('area') != 0:
                err = self.get_error('face_vector_area')
                geo_mean += err[0]
                geo_max = max([geo_max, err[1]])
                err = self.get_error('face_area')
                geo_mean += err[0]
                geo_max = max([geo_max, err[1]])
                n += 2
        if n > 0:
            geo_mean = geo_mean / n
        self.add_error('geometric', geo_mean, geo_mean,
                       self.get_weight('geometric'))

    #--------------------------------------------------------------------------
    #                                   Utilities
    #--------------------------------------------------------------------------

    def axial_forces(self):
        if self.equilibrium != 0:
            return self.mesh.axial_forces()
        else:
            return np.zeros(self.mesh.E)

    def force_resultants(self):
        return self.mesh.force_resultants()

    def applied_loads(self):
        return self.mesh.applied_loads()

    def face_circum_circles(self):
        C = self.face_circumcenters()
        f, v = self.mesh.face_vertices_iterators()
        l = self.mesh.face_lengths()
        r = C[f] - self.mesh.vertices[v]
        d = np.linalg.norm(r, axis=1)
        d = utilities.sum_repeated(d, f)
        r = d / l
        e3 = self.mesh.face_normals()
        e1 = utilities.orthogonal_vectors(e3)
        e2 = np.cross(e3, e1)
        frame = Frame(C, e1, e2, e3)
        circles = Circle(frame, r)
        return circles

    #--------------------------------------------------------------------------
    #                                Errors strings
    #--------------------------------------------------------------------------

    def geometric_error_string(self):
        return self.error_string('geometric')

    def equilibrium_error_string(self):
        return self.error_string('equilibrium')

    def planarity_error_string(self):
        return self.error_string('planarity')


    # -------------------------------------------------------------------------
    #                          Geometric Constraints
    # -------------------------------------------------------------------------

    def normal_constraints(self):
        w = self.get_weight('normal') * self.get_weight('geometric')
        V = self.mesh.V
        F = self.mesh.F
        N = self.N
        f = 3*V + np.arange(F)
        i = np.arange(F)
        i = np.hstack((i, i, i))
        j = np.hstack((f, F+f, 2*F+f))
        X = self.X
        data = 2 * np.hstack((X[f], X[F+f], X[2*F+f])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(F, N))
        r = ((X[f]**2 + X[F+f]**2 + X[2*F+f]**2) + 1) * w
        self.add_iterative_constraint(H, r, 'face_normal_length')

    def planarity_constraints(self):
        w = self.get_weight('planarity')
        V = self.mesh.V
        F = self.mesh.F
        N = self.N
        X = self.X
        f, v1, v2 = self.mesh.face_edge_vertices_iterators(order=True)
        K = f.shape[0]
        f = 3*V + f
        r = ((X[v2] - X[v1]) * X[f] + (X[V+v2] - X[V+v1]) * X[F+f]
             + (X[2*V+v2] - X[2*V+v1]) * X[2*F+f] ) * w
        v1 = np.hstack((v1, V+v1, 2*V+v1))
        v2 = np.hstack((v2, V+v2, 2*V+v2))
        f = np.hstack((f, F+f, 2*F+f))
        i = np.arange(K)
        i = np.hstack((i, i, i, i, i, i, i, i, i))
        j = np.hstack((f, v2, v1))
        data = 2 * np.hstack((X[v2] - X[v1], X[f], -X[f])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(K, N))
        self.add_iterative_constraint(H, r, 'face_normal (planarity)')

    def edge_length_constraints(self):
        w = self.get_weight('edge_length') * self.get_weight('geometric')
        V = self.mesh.V
        E = self.mesh.E
        N = self.N
        N1 = self._N1
        X = self.X
        i = np.arange(E)
        e = N1 + i
        v1, v2 = self.mesh.edge_vertices()
        r = ( X[v1]**2 + X[v2]**2 - 2*X[v1]*X[v2]
            + X[V+v1]**2 + X[V+v2]**2 - 2*X[V+v1]*X[V+v2]
            + X[2*V+v1]**2 + X[2*V+v2]**2 - 2*X[2*V+v1]*X[2*V+v2]
            - X[e]**2 ) * w
        v1 = np.hstack((v1,V+v1,2*V+v1))
        v2 = np.hstack((v2,V+v2,2*V+v2))
        i = np.hstack((i,i,i,i,i,i,i))
        j = np.hstack((v1,v2,e))
        data = 2 * np.hstack((X[v1] - X[v2], X[v2] - X[v1], -X[e])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(E, N))
        self.add_iterative_constraint(H, r, 'edge_length')

    def area_constraints(self):
        w = self.get_weight('area') * self.get_weight('geometric')
        V = self.mesh.V
        F = self.mesh.F
        N = self.N
        N2 = self._N2
        X = self.X
        fi, v1, v2 = self.mesh.face_edge_vertices_iterators(order=True)
        v1a = np.hstack((V+v1,2*V+v1,v1))
        v1b = np.hstack((2*V+v1,v1,V+v1))
        v2a = np.hstack((V+v2,2*V+v2,v2))
        v2b = np.hstack((2*V+v2,v2,V+v2))
        r = 0.5 * (X[v1a] * X[v2b] - X[v1b] * X[v2a]) * w
        k = np.hstack((fi,F+fi,2*F+fi))
        r = utilities.sum_repeated(r,k)
        f = N2 + np.arange(F)
        i = np.hstack((k,k,k,k,np.arange(3*F)))
        j = np.hstack((v1a, v2b, v1b, v2a, f, F+f, 2*F+f))
        data1 = 0.5 * np.hstack((X[v2b], X[v1a], -X[v2a], -X[v1b]))
        data2 = - np.ones(3*F)
        data = np.hstack((data1, data2)) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(3*F, N))
        self.add_iterative_constraint(H, r, 'face_vector_area')

    def vector_area_constraints(self):
        w = self.get_weight('area') * self.get_weight('geometric')
        F = self.mesh.F
        N = self.N
        N2 = self._N2
        f = N2 + np.arange(F)
        i = np.arange(F)
        i = np.hstack((i,i,i,i))
        j = np.hstack((f,F+f,2*F+f,3*F+f))
        X = self.X
        data = 2 * np.hstack((X[f], X[F+f], X[2*F+f], - X[3*F+f])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(F, N))
        r = (X[f]**2 + X[F+f]**2 + X[2*F+f]**2 - X[3*F+f]**2) * w
        self.add_iterative_constraint(H, r, 'face_area')

    def circularity_constraints(self):
        w = self.get_weight('circularity')
        V = self.mesh.V
        F = self.mesh.F
        N = self.N
        N3 = self._N3
        X = self.X
        f, v1, v2 = self.mesh.face_edge_vertices_iterators(order=True)
        K = f.shape[0]
        cx = N3 + np.array(f)
        cy = N3 + np.array(f) + F
        cz = N3 + np.array(f) + 2*F
        v1x = np.array(v1)
        v1y = np.array(V + v1)
        v1z = np.array(2*V + v1)
        v2x = np.array(v2)
        v2y = np.array(V + v2)
        v2z = np.array(2*V + v2)
        jx = np.hstack((v1x, v2x, cx))
        jy = np.hstack((v1y, v2y, cy))
        jz = np.hstack((v1z, v2z, cz))
        datax = np.hstack((X[v1x] - X[cx], -X[v2x] + X[cx], -X[v1x] + X[v2x]))
        rx = 0.5*X[v1x]**2 - 0.5*X[v2x]**2 -X[v1x]*X[cx] + X[v2x]*X[cx]
        datay = np.hstack((X[v1y] - X[cy], -X[v2y] + X[cy], -X[v1y] + X[v2y]))
        ry = 0.5*X[v1y]**2 - 0.5*X[v2y]**2 -X[v1y]*X[cy] + X[v2y]*X[cy]
        dataz = np.hstack((X[v1z] - X[cz], -X[v2z] + X[cz], -X[v1z] + X[v2z]))
        rz = 0.5*X[v1z]**2 - 0.5*X[v2z]**2 -X[v1z]*X[cz] + X[v2z]*X[cz]
        i = np.arange(K)
        i = np.hstack((i,i,i,i,i,i,i,i,i))
        j = np.hstack((jx, jy, jz))
        data = np.hstack((datax, datay, dataz)) * w
        r = (rx + ry + rz) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(K, N))
        self.add_iterative_constraint(H, r,'circularity')

    # -------------------------------------------------------------------------
    #                             Equilibrium
    # -------------------------------------------------------------------------

    def equilibrium_constraints(self):
        w = self.get_weight('equilibrium')
        P = self.mesh.applied_forces
        Fe, Fa = self.mesh.self_weigth_loads
        norm = max(np.max(np.linalg.norm(P, axis=1)), np.max(Fe), np.max(Fa))
        w = 0.1 * w / (norm + 1e-6)
        V = self.mesh.V
        E = self.mesh.E
        F = self.mesh.F
        N = self.N
        N1 = self._N1
        N2 = self._N2
        X = self.X
        constrained = self.mesh.constrained_vertices
        inner = np.invert(np.in1d(np.arange(V), constrained))
        v0, vj = self.mesh.vertex_ring_vertices_iterators(sort=True)
        __, ej = self.mesh.vertex_ring_edges_iterators(sort=True)
        v_mask = np.invert(np.in1d(v0, constrained))
        v0 = v0[v_mask]
        vj = vj[v_mask]
        ej = ej[v_mask]
        ix = np.hstack((v0, v0, v0))
        jx = np.hstack((v0, vj, N1+E+ej))
        datax = np.hstack((X[N1+E+ej], -X[N1+E+ej], X[v0]-X[vj]))
        ri = X[N1+E+ej]*(X[v0]-X[vj])
        ri = utilities.sum_repeated(ri,v0)
        rx = -P[:,0]
        rx[inner] += ri
        iy = ix + V
        jy = np.hstack((V+v0, V+vj, N1+E+ej))
        datay = np.hstack((X[N1+E+ej], -X[N1+E+ej], X[V+v0]-X[V+vj]))
        ri = X[N1+E+ej]*(X[V+v0]-X[V+vj])
        ri = utilities.sum_repeated(ri,v0)
        ry = -P[:,1]
        ry[inner] += ri
        iz = np.hstack((iy+V, 2*V+v0))
        jz = np.hstack((2*V+v0, 2*V+vj, N1+E+ej, N1+ej))
        dataz = np.hstack((X[N1+E+ej], -X[N1+E+ej],
                           X[2*V+v0]-X[2*V+vj], -0.5*Fe[ej]))
        ri = X[N1+E+ej]*(X[2*V+v0]-X[2*V+vj])
        ri = utilities.sum_repeated(ri,v0)
        rz = -P[:,2]
        rz[inner] += ri
        i = np.hstack((ix, iy, iz))
        j = np.hstack((jx, jy, jz))
        data = np.hstack((datax, datay, dataz)) * w
        if self.get_weight('area') != 0:
            vf, fj = self.mesh.vertex_ring_faces_iterators(sort=True)
            f_mask = np.invert(np.in1d(vf, constrained))
            vf = vf[f_mask]
            fj = fj[f_mask]
            Lf = self.mesh.face_lengths()
            iz = 2*V + vf
            jz = N2 + 3*F + fj
            dataz =  - Fa[fj] / Lf[fj] * w
            i = np.hstack((i,iz))
            j = np.hstack((j,jz))
            data = np.hstack((data, dataz))
        r = np.hstack((rx, ry, rz)) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(3*V, N))
        self.add_iterative_constraint(H, r,'equilibrium')

    def boundary_densities_constraints(self):
        w = self.get_weight('equilibrium')
        E = self.mesh.E
        N = self.N
        N1 = self._N1
        X = self.X
        O = N1 + E
        constrained = self.mesh.constrained_vertices
        v1, v2 = self.mesh.edge_vertices()
        mask1 = np.in1d(v1, constrained)
        mask2 = np.in1d(v2, constrained)
        mask = np.logical_and(mask1, mask2)
        #mask = self.mesh.are_boundary_edges()
        j = O + np.arange(E)[mask]
        W = len(j)
        i = np.arange(len(j))
        data = 2*X[j] * w
        r = X[j]**2 * w
        H = sparse.coo_matrix((data,(i,j)), shape=(W,N))
        self.add_constant_constraint(H, r, 'boundary_equilibrium')

    def compression_constraints(self):
        w = self.compression
        sign = np.sign(w)
        E = self.mesh.E
        N = self.N
        N1 = self._N1
        X = self.X
        e = np.arange(E)
        i = np.hstack((e,e))
        j = np.hstack((N1+2*E+e, N1+E+e))
        data = np.hstack((2.0 * X[N1+2*E+e], -sign*np.ones(E))) * abs(w)
        H = sparse.coo_matrix((data,(i,j)), shape=(E,N))
        r = X[N1+2*E+e]**2 * abs(w)
        self.add_iterative_constraint(H, r, 'compression')
        return H, r

    # -------------------------------------------------------------------------
    #                                 Boundary
    # -------------------------------------------------------------------------

    def fixed_boundary_normals_constraints(self):
        w = self.get_weight('fixed_boundary_normals')
        N = self.N
        F = self.mesh.F
        V = self.mesh.V
        f = self.mesh.boundary_faces()
        X = self.X
        W = 3*len(f)
        nx = 3*V + f
        ny = 3*V + F + f
        nz = 3*V + 2*F + f
        data = np.ones(W) * w
        i = np.arange(W)
        j = np.hstack((nx, ny, nz))
        H = sparse.coo_matrix((data,(i,j)), shape=(W,N))
        r = np.hstack((X[nx],X[ny],X[nz])) * w
        self.add_constant_constraint(H, r, 'fixed_boundary_normals')
        return H, r

    # -------------------------------------------------------------------------
    #                                 Build
    # -------------------------------------------------------------------------

    def build_iterative_constraints(self):
        self.add_weight('N', self.N)
        H, r = self.mesh.iterative_constraints(**self.weights)
        self.add_iterative_constraint(H, r, 'mesh_iterative')
        if self.get_weight('planarity') != 0:
            self.normal_constraints()
            self.planarity_constraints()
        if self.get_weight('equilibrium') != 0:
            self.edge_length_constraints()
            self.equilibrium_constraints()
        if self.compression != 0 and self.get_weight('equilibrium') != 0:
            self.compression_constraints()
        if self.get_weight('area') != 0:
            self.area_constraints()
            self.vector_area_constraints()
        if self.get_weight('circularity') != 0:
            self.circularity_constraints()

    def build_constant_constraints(self):
        self.add_weight('N', self.N)
        H, r = self.mesh.constant_constraints(**self.weights)
        self.add_constant_constraint(H, r, 'mesh_constant')
        if self.get_weight('equilibrium') > 0:
            self.boundary_densities_constraints()
        if self.get_weight('fixed_boundary_normals') > 0:
            self.fixed_boundary_normals_constraints()

    def build_constant_fairness(self):
        self.add_weight('N', self.N)
        K, s = self.mesh.fairness_energy(**self.weights)
        self.add_constant_fairness(K, s)

