#!/usr/bin/env python

# -*- coding: utf-8 -*-
__author__ = 'Hui Wang'
#------------------------------------------------------------------------------
import copy

import numpy as np

from scipy import sparse

try:
    from pypardiso import spsolve
except:
    from scipy.sparse.linalg import spsolve

#------------------------------------------------------------------------------
from geometrylab import utilities

#from geometrylab.geometry.meshpy import Mesh # Hui comment

from archgeolab.archgeometry.quadrings import MMesh # Hui

from archgeolab.constraints.constraints_basic import column3D

from archgeolab.constraints.constraints_fairness import con_laplacian_fairness,con_fair_midpoint

#------------------------------------------------------------------------------
""" forcked from geometrylab/optimization/gridshell.py
Hui add: 
    smooth_quadface_data; 
    con_fairness,con_selected_polyline_fairness,con_selected_polylines_fairness,
    con_curve_fold_fairness,con_mesh_fairness,con_boundary_fairness,
    con_tangential_fairness,con_diagonal_mesh_fairness,con_corner_fairness
    
    meshpy.py --> quadrings.py --> gridshell_new.py --> gui_basic.py --> project folder
"""

'''self.mesh calls functions from ProjectFolder: guidedprojection_ + opt_gui_'''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                               Gridshell
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class GridshellNew(MMesh):  # Hui

    def __init__(self, file_name=None):
        #Mesh.__init__(self)
        MMesh.__init__(self)  # Hui

        self._constrained_vertices = None #'boundary' ##Hui comment, otherwise bdry are glided

        self._gliding_vertices = 'constrained'

        self.reference_mesh = None

        self._fixed_vertices = None

        self._beam_load = 0.1 # [N/m]

        self._area_load = 0   # [N/m^2]

        self._vertices_0 = None

        self.boundary_interpolation_N = 8

        self._boundary_polylines = None

        self._gliding_curves = None

        self._applied_forces = None

        self._force_densities = None

        self._handle = np.array([], dtype=int)

        self.__constrained_vertices = np.array([], dtype=int)

        self.__gliding_vertices = np.array([], dtype=int)

        self.__fixed_vertices = np.array([], dtype=int)

        if file_name is not None:
            self.read_obj_file(file_name)

    # -------------------------------------------------------------------------
    #                                 Properties
    # -------------------------------------------------------------------------


    @property
    def beam_load(self):
        return self._beam_load

    @beam_load.setter
    def beam_load(self, load):
        if load != self._beam_load:
            self._beam_load = load
            self.reinitialize_densities()

    @property
    def area_load(self):
        return self._area_load

    @area_load.setter
    def area_load(self, load):
        if load != self._area_load:
            self._area_load = load
            self.reinitialize_densities()

    @property
    def vertices_0(self):
        if self._vertices_0 is None:
            self._vertices_0 = np.copy(self.vertices)
        return np.copy(self._vertices_0)

    @vertices_0.setter
    def vertices_0(self, vertices):
        self._vertices_0 = vertices

    @property
    def handle(self):
        return self._handle

    @handle.setter
    def handle(self, vertex_index):
        if vertex_index is None:
            self._handle = np.array([], dtype=int)
        else:
            self._handle = np.array(vertex_index, dtype=int)

    @property
    def force_densities(self):
        if self._force_densities is None:
            self.lsqr_force_densities()
        densities = np.copy(self._force_densities)
        return densities

    @force_densities.setter
    def force_densities(self, densities):
        self._force_densities = densities

    @property
    def self_weigth_loads(self):
        edge_load, area_load = self._load_array()
        return edge_load, area_load

    @property
    def applied_forces(self):
        if self._applied_forces is None:
            self._applied_forces = np.zeros((self.V, 3), dtype=float)
        forces = np.copy(self._applied_forces)
        return forces

    @applied_forces.setter
    def applied_forces(self, forces):
        self._applied_forces = forces
        self.reinitialize_densities()

    @property
    def axial_forces(self):
        N = self.force_densities * self.edge_lengths()
        return N

    @property
    def fixed_vertices(self):
        if self.__fixed_vertices is None:
            self.__fixed_vertices = np.array([], dtype=int)
        if self._fixed_vertices is None:
            return self.__fixed_vertices
        elif self._fixed_vertices == 'boundary': # Hui: change 'is' to ==
            fix = np.hstack((self.__fixed_vertices, self.boundary_vertices()))
            return np.unique(fix)
        elif self._fixed_vertices == 'constrained':
            fix = np.hstack((self.__fixed_vertices, self.constrained_vertices))
            return np.unique(fix)
        else:
            return self.__fixed_vertices

    @fixed_vertices.setter
    def fixed_vertices(self, vertex_indices):
        if vertex_indices == 'constrained':
            self._fixed_vertices = 'constrained'
        elif vertex_indices == 'boundary':
            self._fixed_vertices = 'boundary'
        elif vertex_indices == None:
            self._fixed_vertices = None
        else:
            self.__fixed_vertices = vertex_indices

    @property
    def constrained_vertices(self):
        if self.__constrained_vertices is None:
            self.__constrained_vertices = np.array([], dtype=int)
        if self._constrained_vertices is None:
            self.__constrained_vertices = np.array([], dtype=int)
            return self.__constrained_vertices
        elif self._constrained_vertices == 'boundary':
            constrained = np.hstack((self.__constrained_vertices,
                                     self.boundary_vertices()))
            return np.unique(constrained)
        else:
            return self.__constrained_vertices

    @constrained_vertices.setter
    def constrained_vertices(self, vertex_indices):
        if vertex_indices == 'boundary':
            self._constrained_vertices = 'boundary'
        elif vertex_indices is None:
            self._constrained_vertices = None
        else:
            self.__constrained_vertices = vertex_indices
        self.reinitialize_densities()

    @property
    def gliding_vertices(self):
        if self.__gliding_vertices is None:
            self.__gliding_vertices = np.array([], dtype=int)
        if self._gliding_vertices is None:
            return self.__gliding_vertices
        elif self._gliding_vertices == 'constrained':
            gli = np.hstack((self.__gliding_vertices, self.constrained_vertices))
            return np.unique(gli)
        elif self._gliding_vertices == 'boundary':
            gli = np.hstack((self.__gliding_vertices, self.boundary_vertices()))
            return np.unique(gli)
        else:
            return self.__gliding_vertices

    @gliding_vertices.setter
    def gliding_vertices(self, gliding_vertices):
        if gliding_vertices == 'constrained':
            self._fixed_vertices = 'constrained'
        elif gliding_vertices == 'boundary':
            self._gliding_vertices = 'boundary'
        elif gliding_vertices is None:
            self._gliding_vertices = None
        else:
            self.__gliding_vertices = gliding_vertices

    @property
    def gliding_curves(self):
        if self._gliding_curves is None:
            self.update_boundary()
        return self._gliding_curves

    # -------------------------------------------------------------------------
    #                                  Copy
    # -------------------------------------------------------------------------

    def copy_gridshell(self):
        copy_gridshell = GridshellNew()
        copy_gridshell.__dict__ = copy.deepcopy(self.__dict__)
        return copy_gridshell

    def import_mesh(self, mesh):
        for key in mesh.__dict__:
            self.__dict__[key] = mesh.__dict__[key]
        self.topology_update()
        self.initialize()

    # -------------------------------------------------------------------------
    #                                  Data
    # -------------------------------------------------------------------------

    def force_resultants(self):
        v, vj = self.vertex_ring_vertices_iterators(sort=True)
        _, ej = self.vertex_ring_edges_iterators(sort=True)
        Vi = self.vertices[v]
        Vj = self.vertices[vj]
        Vij = Vi - Vj
        W = self.force_densities
        Wij = np.array([W[ej], W[ej], W[ej]]).T
        S = Wij*Vij
        S = utilities.sum_repeated(S, v)
        return S

    def loads(self):
        return self._force_array()

    # -------------------------------------------------------------------------
    #                                 Error
    # -------------------------------------------------------------------------

    def equilibrium_error(self):
        F = self._force_array()
        P = self.force_resultants()
        error  = (F + P)
        norm = np.max(np.linalg.norm(F, axis=1))
        if norm == 0:
            norm = 1
        error = np.linalg.norm(error, axis=1) / norm
        error[self.constrained_vertices] = 0
        return error

    # -------------------------------------------------------------------------
    #                                Reading
    # -------------------------------------------------------------------------

    def initialize(self):
        self.set_reference()
        self.set_restart()
        self.set_boundary()

    def reinitialize_densities(self):
        self.force_densities = None

    def topology_update(self):
        super(GridshellNew, self).topology_update()
        self._vertices_0 = None
        self.set_vertices_0()
        self.reset_forces()
        self.reset_fixed()
        self.reset_constraints()
        self.reset_gliding()
        self._gliding_curves = None
        self.reinitialize_densities()

    def set_vertices_0(self):
        self._vertices_0 = np.copy(self.vertices)

    def read_obj_file(self, file_name):
        super(GridshellNew, self).read_obj_file(file_name)
        self.initialize()

    # -------------------------------------------------------------------------
    #                                State
    # -------------------------------------------------------------------------

    def reset(self): # Hui comment: used in windown reset button
        self.vertices = np.copy(self.vertices_0)
        self.set_boundary()
        self.reinitialize_densities()

    def restart(self):
        self.read_obj_file(self.name + '.obj')

    def set_reference(self): # to change: how to add last obj as reference mesh??
        # Hui set:
        # if self.last_object is not None:
        #     self.last_object.geometry.set_reference()
        # else:
        #     self.current_object.geometry.set_reference()
        #print('gridshell')    
        self.reference_mesh = super(GridshellNew, self).copy_mesh()
        self.set_boundary()
        #print(self.reference_mesh.V) # Hui check

    def set_restart(self):
        self.vertices_0 = np.copy(self.vertices)

    # -------------------------------------------------------------------------
    #                                Apply
    # -------------------------------------------------------------------------

    def apply_force(self, force, vertex_indices):
        F = self.applied_forces
        F[vertex_indices,:] += force
        self._applied_forces = F
        self.reinitialize_densities()

    def apply_boundary_load(self, tangent_load, normal_load, vertex_indices):
        F = self.applied_forces
        v = np.array(vertex_indices, dtype=int)
        tangents = self.boundary_tangents(normalize=False)
        A = np.linalg.norm(tangents, axis=1, keepdims=True)
        normals = self.boundary_normals()
        F[v,:] += tangent_load * tangents[v,:]
        F[v,:] += normal_load * A[v] * -normals[v,:]
        self.applied_forces = F
        self.reinitialize_densities()

    def reset_forces(self):
        self._applied_forces = None
        self.reinitialize_densities()

    def constrain(self, vertex_indices):
        v = np.array(vertex_indices, dtype=int)
        self.__constrained_vertices = np.unique(v)
        self.reinitialize_densities()

    def reset_constraints(self):
        self.__constrained_vertices = None
        self.reinitialize_densities()

    def glide(self, vertex_indices):
        v = np.array(vertex_indices, dtype=int)
        self.__gliding_vertices = np.unique(v)

    def reset_gliding(self):
        self.__gliding_vertices = None

    def fix(self, vertex_indices):
        v = np.array(vertex_indices, dtype=int)
        self.__fixed_vertices = np.unique(np.hstack((v, self.__fixed_vertices)))

    def unfix(self, vertex_indices):
        mask = np.in1d(self.__fixed_vertices, vertex_indices)
        mask = np.invert(mask)
        self.__fixed_vertices = self.__fixed_vertices[mask]

    def reset_fixed(self):
        self.__fixed_vertices = None

    def fix_double_boundary(self):
        self.__fixed_vertices = self.double_boundary_vertices()

    # -------------------------------------------------------------------------
    #                                Loads
    # -------------------------------------------------------------------------

    def _load_array(self):
        beam_load = np.ones(self.E) * self.beam_load
        area_load = np.ones(self.F) * self.area_load
        return beam_load, area_load

    def _force_array(self):
        Fe, Fa = self._load_array() #Hui note:  = 0.1, 0
        Le = self.edge_lengths()
        ve, ej = self.vertex_ring_edges_iterators(sort=True)
        vf, fj = self.vertex_ring_faces_iterators(sort=True)
        Fe = - Fe[ej] * Le[ej] / 4 # each edge: -le * 0.1/2, n can==1,2,4
        Fe = utilities.sum_repeated(Fe, ve) #Hui: each vertex: sum edge-Fe
        #print(Fe)
        
        Lf = self.face_lengths() # Fa=0, then Lf no use here
        #print(Lf)
        A = self.face_areas()
        Fa = - Fa[fj] * A[fj] / Lf[fj]
        Fa = utilities.sum_repeated(Fa, vf)
        
        P = np.zeros((self.V, 3))
        P = P + self.applied_forces # Hui note: applied_forces == 0 matrix
        #print(self.applied_forces)
        P[:,2] += (Fe + Fa) # only Z-load, gravity type
        return P

    def lsqr_force_densities(self, dumped=True):
        if self._force_densities is None:
            E = self.E
            constrained = self.constrained_vertices #Hui note: =boundary_vertices()
            v0, vj = self.vertex_ring_vertices_iterators(sort=True)
            v0, ej = self.vertex_ring_edges_iterators(sort=True)
            mask = np.invert(np.in1d(v0, constrained))
            v0 = v0[mask]; vj = vj[mask]; ej = ej[mask]
            V0 = self.vertices[v0,:] #Hui note: =innver v0
            Vj = self.vertices[vj,:] #Hui note: =innver vi
            inner = np.unique(v0)
            W = inner.shape[0]
            iv = utilities.repeated_range(v0) # row index from smallest [0,0,0,0,1...]
            #print(v0,iv,ej)
            i = np.hstack((iv, W + iv, 2*W + iv))
            j = np.hstack((ej, ej, ej))
            data = (V0 - Vj).flatten('F') # note: if Vj-V0, then P=+Fe*le/2, return D no change
            M = sparse.coo_matrix((data,(i,j)), shape=(3*W, E))
            #print(M)
            
            P = self._force_array()[inner] #(Vinner,3)Matrix=[0;0;Fe(at each v)]
            #print(P)
            norm = np.linalg.norm(P, axis=1)
            self._load_scale = 3./(np.max(norm) + 1e-6)
            
            P = P.flatten('F')
            
            if dumped:
                A = (M.T).dot(M) + 1e-6 * sparse.eye(M.shape[1])
                b = (M.T).dot(P)
                D = spsolve(A, b)
            else:
                D = sparse.linalg.lsqr(M, P)[0]
                
            self._force_densities = -D
            #print('fd')

    def densities_space(self, edge_index, closed_boundary=False, fairing=False,
                       curvature_fairing=False):
        E = self.E
        boundary = self.boundary_vertices()
        v0, vj = self.vertex_ring_vertices_iterators(sort=True)
        v0, ej = self.vertex_ring_edges_iterators(sort=True)
        bound = np.in1d(v0, boundary)
        mask = np.invert(bound)
        v0i = v0[mask]; vji = vj[mask]; eji = ej[mask]
        V0 = self.vertices[v0i,:]
        Vj = self.vertices[vji,:]
        inner = np.unique(v0i)
        W = inner.shape[0]
        iv = utilities.repeated_range(v0i)
        i = np.hstack((iv, W + iv, 2*W + iv))
        j = np.hstack((eji, eji, eji))
        data = (V0 - Vj).flatten('F')
        P = np.zeros(3*W)
        M = sparse.coo_matrix((data,(i,j)), shape=(3*W, E))
        if type(edge_index) is not list:
            edge_index = [edge_index]
        W = len(edge_index)
        data = np.ones(W)
        i = np.zeros(W)
        j = np.array(edge_index)
        Pi = np.array([W])
        Mi = sparse.coo_matrix((data,(i,j)), shape=(1, E))
        M = sparse.vstack((M, Mi))
        P = np.hstack((P, Pi))
        if closed_boundary:
            curves = self.boundary_curves(corner_split=False)
            for curve in curves:
                mask = np.invert(np.in1d(v0, curve))
                v0c = v0[mask]; vjc = vj[mask]; ejc = ej[mask]
                V0 = self.vertices[v0c,:]
                Vj = self.vertices[vjc,:]
                W = len(v0c)
                i = np.zeros(W)
                i = np.hstack((i, i+1, i+2))
                j = np.hstack((ejc, ejc, ejc))
                data = (V0 - Vj).flatten('F') * 10
                Mc = sparse.coo_matrix((data,(i,j)), shape=(3, E))
                Pc = np.zeros(3)
                M = sparse.vstack((M, Mc))
                P = np.hstack((P, Pc))
        if fairing or curvature_fairing:
            e, e1, e2, e11, e22 = self.edge_side_edges()
            W = len(e)
            i = np.arange(W)
            w = 0.02
            one = np.ones(W)
            i = np.hstack((i, i, i, i, i))
            j = np.hstack((e, e1, e2, e11, e22))
            Pi = np.zeros(W)
        if fairing:
            data = np.hstack((3*one, -2*one, -2*one, .5*one, .5*one)) * w
            Mi = sparse.coo_matrix((data,(i,j)), shape=(W, E))
            M = sparse.vstack((M, Mi))
            P = np.hstack((P, Pi))
        if curvature_fairing:
            data = np.hstack((-7*one, 3*one, 3*one, .5*one, .5*one)) * w
            Mi = sparse.coo_matrix((data,(i,j)), shape=(W, E))
            M = sparse.vstack((M, Mi))
            P = np.hstack((P, Pi))
        A = (M.T).dot(M) + 1e-10 * sparse.eye(M.shape[1])
        b = (M.T).dot(P)
        D = spsolve(A, b)
        #D = sparse.linalg.lsqr(M, P)[0]
        self._force_densities = -D
        self.extract_boundary_load()

    def extract_boundary_load(self):
        F = self.force_resultants()
        b = self.boundary_vertices()
        self.applied_forces[b,:] = F[b,:]

    # -------------------------------------------------------------------------
    #                               Boundary
    # -------------------------------------------------------------------------

    def _interpolate_boundary(self):
        N = self.boundary_interpolation_N
        polylines = self.boundary_polylines()
        if N > 0:
            for polyline in polylines:
                polyline.refine(N)
        self._boundary_polylines = polylines

    def set_boundary(self):
        self._interpolate_boundary()
        self.update_boundary()

    def update_boundary(self):
        curves = self.boundary_curves(corner_split=False)
        gliding_curves = []
        if self.gliding_vertices.shape[0] > 0:
            for curve in curves:
                gliding = np.in1d(curve, self.gliding_vertices)
                gliding_curve = curve[gliding]
                gliding_curves.append(np.unique(gliding_curve))
        self._gliding_curves = gliding_curves

    # -------------------------------------------------------------------------
    #                             Fixed and Gliding
    # -------------------------------------------------------------------------

    def fixed_vertices_constraints(self, **kwargs):
        N = kwargs.get('N', 3*self.V)
        if 'fixed_vertices' in kwargs:
            w = kwargs.get('fixed_vertices')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        fixed_vertices = np.unique(np.hstack((self.fixed_vertices, self.handle)))
        W = len(fixed_vertices)
        v = np.array(fixed_vertices)
        j = np.hstack((v,V+v,2*V+v))
        i = np.arange(W)
        i = np.hstack((i,W+i,2*W+i))
        r = self.vertices[v,:]
        r = np.reshape(r, 3*W, order='F') * w
        data = np.ones([3*W]) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(3*W, N))
        if 'settings' in kwargs:
            settings = kwargs['settings']
            settings.add_constraint('fixed_vertices', 3*W)
            #print(settings)
        #print('fix_weight=',w) # default setting w=10
        return H, r

    def self_closeness_constraints(self, **kwargs):
        N = kwargs.get('N', 3*self.V)
        if 'self_closeness' in kwargs:
            w = kwargs.get('self_closeness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        #fixed_vertices = np.arange(self.V)
        W = self.V#len(fixed_vertices) # Hui simplify
        v = np.arange(W)# np.array(fixed_vertices)
        j = np.hstack((v,V+v,2*V+v))
        i = np.arange(W)
        i = np.hstack((i,W+i,2*W+i))
        ver = self.vertices#[v,:]
        r = np.reshape(ver, 3*W, order='F') * w
        data = np.ones([3*W]) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(3*W, N))
        if 'settings' in kwargs:
            settings = kwargs['settings']
            settings.add_constraint('self_closeness', 3*W)
        return H, r

    def fixed_corners_constraints(self, **kwargs):
        N = kwargs.get('N', 3*self.V)
        if 'fixed_corners' in kwargs:
            w = kwargs.get('fixed_corners')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        fixed_vertices = self.mesh_corners()
        W = len(fixed_vertices)
        v = np.array(fixed_vertices)
        j = np.hstack((v,V+v,2*V+v))
        i = np.arange(W)
        i = np.hstack((i,W+i,2*W+i))
        r = np.reshape(self.vertices[v], 3*W, order='F') * w
        data = np.ones(3*W) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(3*W, N))
        if 'settings' in kwargs:
            settings = kwargs['settings']
            settings.add_constraint('fixed_corners', 3*W)
        return H, r

    def gliding_constraints(self, **kwargs):
        N = kwargs.get('N', 3*self.V)
        if 'gliding' in kwargs:
            w = kwargs.get('gliding')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        gliding_curves = self.gliding_curves
        H = sparse.coo_matrix(([0],([0],[0])), shape=(1, N))
        r = np.array([0])
        for b in range(len(gliding_curves)):
            curve = gliding_curves[b]
            W = len(curve)
            polyline = self._boundary_polylines[b]
            closest = polyline.closest_vertices(self.vertices[curve,:])
            points = polyline.vertices[closest,:]
            t = polyline.vertex_tangents()[closest]
            n1 = utilities.orthogonal_vectors(t)
            n2 = np.cross(t,n1)
            k = np.arange(W)
            r1 = w * np.einsum('ij,ij->i', points, n1)
            data1 = w * n1.flatten('F')
            r2 = w * np.einsum('ij,ij->i', points, n2)
            data2 = w * n2.flatten('F')
            i = np.hstack((k,k,k,W+k,W+k,W+k))
            j = np.hstack((curve,V+curve,2*V+curve,curve,V+curve,2*V+curve))
            data = np.hstack((data1,data2))
            r = np.hstack((r,r1,r2))
            Hb = sparse.coo_matrix((data,(i,j)), shape=(2*W, N))
            H = sparse.vstack([H,Hb])
            if 'settings' in kwargs:
                settings = kwargs['settings']
                settings.add_constraint('gliding_boundary_' + str(b), 2*W)
        return H, r

    def boundary_closeness_constraints(self, **kwargs):
        N = kwargs.get('N', 3*self.V)
        if 'boundary_closeness' in kwargs:
            w = kwargs.get('boundary_closeness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        gliding_curves = self.boundary_curves(corner_split=True)
        H = sparse.coo_matrix(([0],([0],[0])), shape=(1, N))
        r = np.array([0])
        for b in range(len(gliding_curves)):
            curve = gliding_curves[b]
            W = len(curve)
            polyline = self._boundary_polylines[b]
            closest = polyline.closest_vertices(self.vertices[curve,:])
            points = polyline.vertices[closest,:]
            t = polyline.vertex_tangents()[closest]
            n1 = utilities.orthogonal_vectors(t)
            n2 = np.cross(t,n1)
            k = np.arange(W)
            r1 = w * np.einsum('ij,ij->i', points, n1)
            data1 = w * n1.flatten('F')
            r2 = w * np.einsum('ij,ij->i', points, n2)
            data2 = w * n2.flatten('F')
            i = np.hstack((k, k, k, W+k, W+k, W+k))
            j = np.hstack((curve,V+curve,2*V+curve,curve,V+curve,2*V+curve))
            data = np.hstack((data1, data2))
            r = np.hstack((r, r1, r2))
            Hb = sparse.coo_matrix((data,(i,j)), shape=(2*W, N))
            H = sparse.vstack([H, Hb])
            if 'settings' in kwargs:
                settings = kwargs['settings']
                settings.add_constraint('closeness_boundary_' + str(b), 2*W)
        return H, r

    # -------------------------------------------------------------------------
    #                                  Reference
    # -------------------------------------------------------------------------

    def reference_closeness_constraints(self, **kwargs):
        """
        Hui note: 
            (vi-pi) * ni = 0
        vi: variables
        pi: computed closest point
        ni: computed normal
        """
        N = kwargs.get('N', 3*self.V)
        if 'reference_closeness' in kwargs:
            w = kwargs.get('reference_closeness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        reference_mesh = self.reference_mesh
        #print('con_ref',V,self.reference_mesh.V) # Hui check
        normals = reference_mesh.vertex_normals()
        closest = reference_mesh.closest_vertices(self.vertices)
        points = reference_mesh.vertices[closest,:]
        normals = normals[closest,:]
        r = w * np.einsum('ij,ij->i', points, normals)
        data = w * normals.flatten('F')
        v = np.arange(V)
        j = np.hstack((v, V+v, 2*V+v))
        i = np.hstack((v, v, v))
        H = sparse.coo_matrix((data, (i,j)), shape=(V, N))
        if 'settings' in kwargs:
            settings = kwargs['settings']
            settings.add_constraint('reference_closeness', V)
        return H, r

    # -------------------------------------------------------------------------
    #                                 Fairness
    # -------------------------------------------------------------------------

    def mesh_fairness(self, **kwargs):
        """ Hui note: 
            if n=4: 
                v1+v3 = 2*v; v2+v4 = 2*v
            else:
                v1+...vn = n*v
        data1 = -1/n
        data2 = 1
        data3 = data4 = -0.5
        data5 = data6 = 1
        """
        N = kwargs.get('N', 3*self.V)
        if 'mesh_fairness' in kwargs:
            w = kwargs.get('mesh_fairness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        bound = self.boundary_vertices()
        v0, vj, l = self.vertex_ring_vertices_iterators(order=True,
                                                     return_lengths=True)
        inner = np.invert(np.in1d(v0, bound))
        v0 = v0[inner]
        vj = vj[inner]
        quad = np.in1d(v0, np.where(l == 4)[0])
        other = np.invert(quad)
        v0o = v0[other]
        if v0o.shape[0] > 0:
            i1 = utilities.repeated_range(v0o)
            j1 = vj[other]
            data1 = -np.ones(v0o.shape[0]) / l[v0o]
            off1 = np.amax(i1) + 1
            i2 = np.arange(off1)
            j2 = np.unique(v0o)
            data2 = np.ones(j2.shape[0])
        else:
            i1 = i2 = j1= j2 = data1 = data2 = np.array([])
            off1 = 0
        v0q = v0[quad][0::2]
        if v0q.shape[0] > 0:
            i3 = utilities.repeated_range(v0q, offset=off1)
            j3 = vj[quad][0::2]
            data3 = np.repeat(-.5,v0q.shape[0]) # =[-0.5,-0.5,...]
            i4 = i3 + np.amax(i3) + 1 - off1
            j4 = vj[quad][1::2]
            data4 = data3
            i5 = np.arange(off1, off1 + v0q.shape[0])
            j5 = np.unique(v0q)
            data5 = np.ones(i5.shape[0])
        else:
            i3 = i4 = i5 = j3 = j4 = j5 = data3 = data4 = data5 = np.array([])
        i = np.hstack((i1,i2,i3,i4,i5))
        j = np.hstack((j1,j2,j3,j4,j5,j5))
        data = w * np.hstack((data1, data2, data3, data4, data5))
        W = int(np.amax(i))
        Kx = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        j = V+j
        Ky = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        j = V+j
        Kz = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        K = sparse.vstack([Kx, Ky, Kz])
        s = np.zeros(int(3*(W+1)))
        #np.savetxt('test.txt', K.toarray(),fmt='% 2d')
        return K, s

    def tangential_fairness(self, **kwargs):
        """ Hui note same as mesh_fairness but data * ti: 
            if n=4: 
                v1+v3 = 2*v; v2+v4 = 2*v
            else:
                v1+...vn = n*v
        data1 = -1/n
        data2 = 1
        data3 = data4 = -0.5
        data5 = data6 = 1      
          
        above matrix K
        K1 = K * t1
        K2 = K * t2
        """
        N = kwargs.get('N', 3*self.V)
        if 'tangential_fairness' in kwargs:
            w = kwargs.get('tangential_fairness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        normals = self.vertex_normals()
        t1 = utilities.orthogonal_vectors(normals)
        t2 = np.cross(normals,t1)
        inner = self.inner_vertices()
        v0, vj, l = self.vertex_ring_vertices_iterators(order=True,
                                                     return_lengths=True)
        inner = np.in1d(v0, inner)
        v0 = v0[inner]
        vj = vj[inner]
        quad = np.in1d(v0, np.where(l == 4)[0])
        other = np.invert(quad)
        v0o = v0[other]
        if v0o.shape[0] > 0:
            i1 = utilities.repeated_range(v0o)
            j1 = vj[other]
            data1 = -np.ones(v0o.shape[0]) / l[v0o]
            off1 = np.amax(i1) + 1
            i2 = np.arange(off1)
            j2 = np.unique(v0o)
            data2 = np.ones(j2.shape[0])
        else:
            i1 = i2 = j1= j2 = data1 = data2 = []
            off1 = 0
        v0q = v0[quad][0::2]
        if v0q.shape[0] > 0:
            i3 = utilities.repeated_range(v0q, offset=off1)
            j3 = vj[quad][0::2]
            data3 = np.repeat(-.5,v0q.shape[0])
            i4 = i3 + np.amax(i3) + 1 - off1
            j4 = vj[quad][1::2]
            data4 = data3
            i5 = np.arange(off1, off1 + v0q.shape[0])
            j5 = np.unique(v0q)
            data5 = np.ones(i5.shape[0])
        else:
            i3 = i4 = i5 = j3 = j4 = j5 = data3 = data4 = data5 = []
        i = np.hstack((i1,i2,i3,i4,i5))
        j = np.hstack((j1,j2,j3,j4,j5,j5))
        data = w * np.hstack((data1, data2, data3, data4, data5))
        t1=np.hstack((t1[v0o,0],t1[j2,0],t1[v0q,0],t1[v0q,0],t1[j5,0],t1[j5,0],
                     t1[v0o,1],t1[j2,1],t1[v0q,1],t1[v0q,1],t1[j5,1],t1[j5,1],
                     t1[v0o,2],t1[j2,2],t1[v0q,2],t1[v0q,2],t1[j5,2],t1[j5,2]))
        t2=np.hstack((t2[v0o,0],t2[j2,0],t2[v0q,0],t2[v0q,0],t2[j5,0],t2[j5,0],
                     t2[v0o,1],t2[j2,1],t2[v0q,1],t2[v0q,1],t2[j5,1],t2[j5,1],
                     t2[v0o,2],t2[j2,2],t2[v0q,2],t2[v0q,2],t2[j5,2],t2[j5,2]))
        W = int(np.amax(i))
        i = np.hstack((i,i,i))
        j = np.hstack((j,V+j,2*V+j))
        data1 = np.hstack((data,data,data)) * t1
        K1 = sparse.coo_matrix((data1,(i,j)), shape=(W+1,N))
        data2 = np.hstack((data,data,data)) * t2
        K2 = sparse.coo_matrix((data2,(i,j)), shape=(W+1,N))
        K = sparse.vstack([K1,K2])
        s = np.zeros(int(2*(W+1)))
        return K, s

    def boundary_fairness(self, **kwargs):
        "Hui note: Vl+Vr = 2*V "
        N = kwargs.get('N', 3*self.V)
        if 'boundary_fairness' in kwargs:
            w = kwargs.get('boundary_fairness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        boundaries = self.boundary_curves(corner_split=True)
        corners = self.mesh_corners()
        K = sparse.coo_matrix(([0],([0],[0])), shape=(1,N))
        for boundary in boundaries:
            mask = np.invert(np.in1d(boundary, corners))
            v0 = boundary[mask]
            W = v0.shape[0]
            vn = np.roll(boundary,1)[mask]
            vp = np.roll(boundary,-1)[mask]
            one = np.ones(W)
            v = np.arange(W)
            i = np.hstack((v,v,v,W+v,W+v,W+v,2*W+v,2*W+v,2*W+v))
            j = np.hstack((v0,vp,vn,V+v0,V+vp,V+vn,2*V+v0,2*V+vp,2*V+vn))
            data = np.hstack((2*one,-one,-one))
            data = np.hstack((data,data,data)) * w
            Ki = sparse.coo_matrix((data,(i,j)), shape=(3*W,N))
            K = sparse.vstack([K, Ki])
        s = np.zeros((K.shape[0]))
        return K, s

    def spring_fairness(self, **kwargs):
        "Hui note: edges(V1V2): V1~=V2 "
        N = kwargs.get('N', 3*self.V)
        if 'spring_fairness' in kwargs:
            w = kwargs.get('spring_fairness')
        else:
            w = kwargs.get('w', 1)
        V = self.V
        E = self.E
        v1, v2 = self.edge_vertices()
        i = np.arange(E)
        i = np.hstack((i,i,E+i,E+i,2*E+i,2*E+i))
        j = np.hstack((v1,v2,V+v1,V+v2,2*V+v1,2*V+v2))
        one = np.ones(E)
        data = np.hstack((one, -one, one, -one, one, -one)) * w
        K = sparse.coo_matrix((data,(i,j)), shape=(3*E,N))
        s = np.zeros((3*E))
        return K, s


    # -------------------------------------------------------
    #                   Hui add fairness functions
    # -------------------------------------------------------   
    def con_fairness(self,move,N,poly_edges,
                     eboundary3=2,einner3=2,einner5=2,einner6=2,ecorner4=0,
                     **kwargs):
        """current: for valence 4 star: inner non-sigular points
            if self.vertices, move = 0
            elif velocity vectors, move = numvlc
            elif refer mesh, move = Nr-3*V
            2*v-v1-v3=0; 2*v-v2-v4=0
        """
        N = kwargs.get('N')
        w = kwargs.get('fairness')
        Vnum = self.V
        v,v1,v2,v3,v4 = self.ver_regular_star.T

        if kwargs.get('is_fair_curve_fold'): #see con_selected_polyline_fairness / con_curve_fold_fairness.
            # ii = self.get_polylines_fair_index(self.poly_edges,vvv,vva,vvb)
            # print('--',self.poly_edges,vvv[ii])
            # jj = np.setdiff1d(np.arange(len(vvv)), ii)
            # vvv,vva,vvb = vvv[jj],vva[jj],vvb[jj]
            poly,_,_ = self.get_polylines_from_edges(poly_edges)
            ii = np.invert(np.in1d(v,poly))
            v,v1,v2,v3,v4 = v[ii], v1[ii], v2[ii], v3[ii], v4[ii]            
        vvv = np.r_[v, v]
        vva = np.r_[v1,v2]
        vvb = np.r_[v3,v4]          
        num = len(vvv)      
        arr = np.arange(3*num)
        one = np.ones(3*num)
        row = np.tile(arr,3)
        
        cv = column3D(vvv, move,Vnum)
        ca = column3D(vva,move,Vnum)
        cb = column3D(vvb,move,Vnum)
        col = np.r_[cv,ca,cb]
        data = np.r_[2*one,-one,-one]
        K = sparse.coo_matrix((data,(row,col)), shape=(3*num, N))

        if eboundary3:
            "vl+vr = 2*v"
            try:
                v,vl,vr = self.ver_bod_valence3_neib
                
                if kwargs.get('is_fair_curve_fold'):
                    ii = np.invert(np.in1d(v,poly))
                    v,vl,vr = v[ii], vl[ii], vr[ii]
                
                K3 = con_fair_midpoint(v,vl,vr,move,Vnum,N)*eboundary3
                K = sparse.vstack((K,K3)) 
                v,vl,vr = self.vertex_corner_valence3_neib()
                Kc = con_fair_midpoint(v,vl,vr,move,Vnum,N)*eboundary3*1
                K = sparse.vstack((K,Kc)) 
            except:
                pass
        if einner3:
            "v1+v2+v3 = 3*v"
            try:
                v,neib = self.vertex_inn_valence3_neib
                K3 = con_laplacian_fairness(v,neib,move,Vnum,N)*einner3  
                K = sparse.vstack((K,K3)) 
            except:
                pass
        if einner5:
            "v1+v2+v3+v4+v5 = 5*v"
            try:
                v,neib = self.vertex_inn_valence5_neib
                K5 = con_laplacian_fairness(v,neib,move,Vnum,N)*einner5  
                K = sparse.vstack((K,K5)) 
            except:
                pass
        if einner6:
            "v1+v2+v3+v4+v5+v6 = 6*v"
            try:
                v,neib = self.vertex_inn_valence6_neib
                K6 = con_laplacian_fairness(v,neib,move,Vnum,N)*einner6
                K = sparse.vstack((K,K6)) 
            except:
                pass
        if ecorner4:
            "v1+---+v8 = 8*v"
            try:
                v,neib = self.ver_corner_valence4_neib
                Kc4 = con_laplacian_fairness(v,neib,move,Vnum,N)*ecorner4
                K = sparse.vstack((K,Kc4)) 
            except:
                pass         
        s = np.zeros(K.shape[0]) 
        return K  * w, s
    
    def con_selected_polyline_fairness(self,move,N,poly_edges,**kwargs):
        N = kwargs.get('N')
        w = kwargs.get('poly_fairness')
        Vnum = self.V
        v,v1,v2,v3,v4 = self.ver_regular_star.T
        vvv = np.r_[v, v]
        vva = np.r_[v1,v2]
        vvb = np.r_[v3,v4]           

        ii = self.get_polylines_fair_index(poly_edges,vvv,vva,vvb)
        vvv,vva,vvb = vvv[ii],vva[ii],vvb[ii]

        num = len(vvv)      
        arr = np.arange(3*num)
        one = np.ones(3*num)
        row = np.tile(arr,3)
        
        cv = column3D(vvv,move,Vnum)
        ca = column3D(vva,move,Vnum)
        cb = column3D(vvb,move,Vnum)
        col = np.r_[cv,ca,cb]
        data = np.r_[2*one,-one,-one]
        K = sparse.coo_matrix((data,(row,col)), shape=(3*num, N))  
        s = np.zeros(K.shape[0]) 
        return K  * w, s
    
    def con_selected_polylines_fairness(self,move,N,poly_edges,**kwargs):
        N = kwargs.get('N')
        w = kwargs.get('polys_fairness')
        def get_inside_polyline_neib(poly):
            v,v1,v2,v3,v4 = self.ver_regular_star.T
            vb,vbl,vbr = self.ver_bod_valence3_neib   
            vv,vvl,vvr = np.r_[v,v,vb], np.r_[v1,v2,vbl], np.r_[v3,v4,vbr]
            vp,vpl,vpr = [],[],[]
            for i in range(len(poly)):
                k = np.where(vv==poly[i])[0]
                for j in k:
                    if vvl[j] in poly and vvr[j] in poly:
                        vp.append(vv[j])
                        vpl.append(vvl[j])
                        vpr.append(vvr[j])
            return np.array(vp),np.array(vpl),np.array(vpr)
        vp,vpl,vpr = get_inside_polyline_neib(poly_edges)
        Vnum = self.V
        num = len(vp)      
        arr = np.arange(3*num)
        one = np.ones(3*num)
        row = np.tile(arr,3)
        
        cv = column3D(vp,move,Vnum)
        ca = column3D(vpl,move,Vnum)
        cb = column3D(vpr,move,Vnum)
        col = np.r_[cv,ca,cb]
        data = np.r_[2*one,-one,-one]
        K = sparse.coo_matrix((data,(row,col)), shape=(3*num, N))  
        s = np.zeros(K.shape[0]) 
        return K  * w, s
    
    def con_curve_fold_fairness(self,poly_edges,**kwargs):
        N = kwargs.get('N')
        if poly_edges is not None:
            poly,_,_ = self.get_polylines_from_edges(poly_edges)
        else:
            poly = None
        K = sparse.coo_matrix(([0],([0],[0])), shape=(1,N))
        s = np.zeros(1)
        if kwargs.get('fair_inner'):
            K1,s1 = self.con_mesh_fairness(poly=poly,**kwargs)
            K = sparse.vstack((K,K1)) 
            s = np.r_[s,s1]
        if kwargs.get('fair_boundary'):
            K2,s2 = self.con_boundary_fairness(poly=poly,**kwargs)
            K = sparse.vstack((K,K2)) 
            s = np.r_[s,s2]
        if kwargs.get('fair_tangential'):
            K3,s3 = self.con_tangential_fairness(poly=poly,**kwargs)
            K = sparse.vstack((K,K3)) 
            s = np.r_[s,s3]
        return K,s
  
    def con_mesh_fairness(self,move=0,poly=None,poly_fair=True,**kwargs):
        """
        refer: gridshell/mesh_fairness
        remove poly vertices
        """
        N = kwargs.get('N')
        w = kwargs.get('mesh_fairness')
        V = self.V
        v0, vj, l = self.vertex_ring_vertices_iterators(order=True,
                                                     return_lengths=True)
        
        if poly is not None:
            ind = np.invert(np.in1d(v0,poly))
            v0 = v0[ind]
            vj = vj[ind]
        
        if False:
            "not includ boundary fairness of valence 4"
            bound = self.boundary_vertices()
            inner = np.invert(np.in1d(v0, bound))
            v0 = v0[inner]
            vj = vj[inner]
        
        quad = np.in1d(v0, np.where(l == 4)[0])
        other = np.invert(quad)
        v0o,vjo = v0[other],vj[other]
        
        if v0o.shape[0] > 0:
            i1 = utilities.repeated_range(v0o)
            j1 = vjo
            data1 = -np.ones(v0o.shape[0]) / l[v0o]
            off1 = np.amax(i1) + 1
            i2 = np.arange(off1)
            j2 = np.unique(v0o)
            data2 = np.ones(j2.shape[0])
        else:
            i1 = i2 = j1= j2 = data1 = data2 = np.array([])
            off1 = 0
        v0q = v0[quad][0::2]
        if v0q.shape[0] > 0:
            i3 = utilities.repeated_range(v0q, offset=off1)
            j3 = vj[quad][0::2]
            data3 = np.repeat(-.5,v0q.shape[0]) # =[-0.5,-0.5,...]
            i4 = i3 + np.amax(i3) + 1 - off1
            j4 = vj[quad][1::2]
            data4 = data3
            i5 = np.arange(off1, off1 + v0q.shape[0])
            j5 = np.unique(v0q)
            data5 = np.ones(i5.shape[0])
        else:
            i3 = i4 = i5 = j3 = j4 = j5 = data3 = data4 = data5 = np.array([])
        i = np.hstack((i1,i2,i3,i4,i5)) 
        j = np.hstack((j1,j2,j3,j4,j5,j5)) + move
        data = np.hstack((data1, data2, data3, data4, data5)) * w
        W = int(np.amax(i))
        Kx = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        j = V+j
        Ky = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        j = V+j
        Kz = sparse.coo_matrix((data,(i,j)), shape=(W+1,N))
        K = sparse.vstack([Kx, Ky, Kz])
        s = np.zeros(int(3*(W+1)))
        
        if poly is not None and poly_fair:
            #Kp,sp = self.con_selected_polyline_fairness(0,N,efair*1) # need check weight
            Kp,sp = self.con_selected_polylines_fairness(0,N,poly,w*1)
            K = sparse.vstack((K,Kp)) 
            s = np.r_[s,sp]
        return K, s   
    
    def con_boundary_fairness(self, move=0,poly=None,**kwargs):
        """
        refer: gridshell/boundary_fairness
        remove poly vertices
        """
        N = kwargs.get('N')
        w = kwargs.get('boundary_fairness')
        V = self.V
        boundaries = self.boundary_curves(corner_split=True)
        corners = self.mesh_corners()
        K = sparse.coo_matrix(([0],([0],[0])), shape=(1,N))
        for boundary in boundaries:
            
            if poly is not None:
                #ind = np.invert(np.in1d(boundary,poly))
                #boundary = boundary[ind]       
                corners = np.r_[corners,poly]
            
            mask = np.invert(np.in1d(boundary, corners))
            v0 = boundary[mask]
            W = v0.shape[0]
            vn = np.roll(boundary,1)[mask]
            vp = np.roll(boundary,-1)[mask]
            one = np.ones(W)
            v = np.arange(W)
            i = np.hstack((v,v,v,W+v,W+v,W+v,2*W+v,2*W+v,2*W+v))
            j = np.hstack((v0,vp,vn,V+v0,V+vp,V+vn,2*V+v0,2*V+vp,2*V+vn)) + move
            data = np.hstack((2*one,-one,-one))
            data = np.hstack((data,data,data))
            Ki = sparse.coo_matrix((data,(i,j)), shape=(3*W,N))
            K = sparse.vstack([K, Ki])
        s = np.zeros((K.shape[0]))
        return K*w, s 

    def con_tangential_fairness(self,move=0,poly=None,**kwargs):
        """
        refer: gridshell/tangential_fairness
        remove poly vertices
        """
        N = kwargs.get('N')
        w = kwargs.get('tangential_fairness')
        V = self.V
        normals = self.vertex_normals()
        t1 = utilities.orthogonal_vectors(normals)
        t2 = np.cross(normals,t1)
        inner = self.inner_vertices()
        v0, vj, l = self.vertex_ring_vertices_iterators(order=True,
                                                     return_lengths=True)
        
        if poly is not None:
            ind = np.invert(np.in1d(v0,poly))
            v0 = v0[ind]
            vj = vj[ind]         
        
        inner = np.in1d(v0, inner)
        v0 = v0[inner]
        vj = vj[inner]
        quad = np.in1d(v0, np.where(l == 4)[0])
        other = np.invert(quad)
        v0o, vjo = v0[other], vj[other]
        if v0o.shape[0] > 0:
            i1 = utilities.repeated_range(v0o)
            j1 = vjo
            data1 = -np.ones(v0o.shape[0]) / l[v0o]
            off1 = np.amax(i1) + 1
            i2 = np.arange(off1)
            j2 = np.unique(v0o)
            data2 = np.ones(j2.shape[0])
        else:
            i1 = i2 = j1= j2 = data1 = data2 = []
            off1 = 0
        v0q = v0[quad][0::2]
        if v0q.shape[0] > 0:
            i3 = utilities.repeated_range(v0q, offset=off1)
            j3 = vj[quad][0::2]
            data3 = np.repeat(-.5,v0q.shape[0])
            i4 = i3 + np.amax(i3) + 1 - off1
            j4 = vj[quad][1::2]
            data4 = data3
            i5 = np.arange(off1, off1 + v0q.shape[0])
            j5 = np.unique(v0q)
            data5 = np.ones(i5.shape[0])
        else:
            i3 = i4 = i5 = j3 = j4 = j5 = data3 = data4 = data5 = []
        i = np.hstack((i1,i2,i3,i4,i5))
        j = np.hstack((j1,j2,j3,j4,j5,j5))
        data = np.hstack((data1, data2, data3, data4, data5))
        t1=np.hstack((t1[v0o,0],t1[j2,0],t1[v0q,0],t1[v0q,0],t1[j5,0],t1[j5,0],
                     t1[v0o,1],t1[j2,1],t1[v0q,1],t1[v0q,1],t1[j5,1],t1[j5,1],
                     t1[v0o,2],t1[j2,2],t1[v0q,2],t1[v0q,2],t1[j5,2],t1[j5,2]))
        t2=np.hstack((t2[v0o,0],t2[j2,0],t2[v0q,0],t2[v0q,0],t2[j5,0],t2[j5,0],
                     t2[v0o,1],t2[j2,1],t2[v0q,1],t2[v0q,1],t2[j5,1],t2[j5,1],
                     t2[v0o,2],t2[j2,2],t2[v0q,2],t2[v0q,2],t2[j5,2],t2[j5,2]))
        W = int(np.amax(i))
        i = np.hstack((i,i,i))
        j = np.hstack((j,V+j,2*V+j)) + move
        data1 = np.hstack((data,data,data)) * t1
        K1 = sparse.coo_matrix((data1,(i,j)), shape=(W+1,N))
        data2 = np.hstack((data,data,data)) * t2
        K2 = sparse.coo_matrix((data2,(i,j)), shape=(W+1,N))
        K = sparse.vstack([K1,K2])
        s = np.zeros(int(2*(W+1)))
        return K*w, s

    def con_diagonal_mesh_fairness(self, **kwargs):
        """mesh_fairness for diagonal mesh
           a   1    d
           2   v    4
           b   3    c
        fairness:
            a+c-2v=0; b+d-2v=0.
        """
        N = kwargs.get('N')
        w = kwargs.get('fairness_diagmesh')
        v,va,vb,vc,vd = self.rr_star_corner
        K13 = con_fair_midpoint(v,va,vc,0,self.V,N)
        K24 = con_fair_midpoint(v,vb,vd,0,self.V,N)
        K = sparse.vstack([K13,K24])
        s = np.zeros(K.shape[0])
        return K*w, s

    def con_corner_fairness(self,move=0,**kwargs):
        N = kwargs.get('N')
        w = kwargs.get('corner_fairness')
        Vnum = self.V
        vl,vc,vr = self.ver_corner_neib
        K = con_fair_midpoint(vc,vl,vr,move,Vnum,N)
        try:
            "v1+---+v8 = 8*v"
            v,neib = self.ver_corner_valence4_neib
            Kc4 = con_laplacian_fairness(v,neib,move,Vnum,N)
            K = sparse.vstack((K,Kc4)) 
        except:
            pass         
        s = np.zeros(K.shape[0]) 
        #self.add_iterative_constraint(K*w, s, 'fairness_diagmesh')    
        return K*w,s

    
    # -------------------------------------------------------------------------
    #                                 Build
    # -------------------------------------------------------------------------

    def constant_constraints(self, **kwargs):
        N = kwargs.get('N')
        i = j = data = r = np.zeros([0])
        H = sparse.coo_matrix((data,(i,j)), shape=(0,N))
        if kwargs.get('fixed_vertices') != 0 and kwargs.get('fixed_vertices') is not None:
            if len(self.fixed_vertices) + len(self.handle) != 0:
                Hi, ri = self.fixed_vertices_constraints(**kwargs)
                H = sparse.vstack([H, Hi])
                r = np.hstack((r, ri))
                #print('1')
        if kwargs.get('fixed_corners') != 0 and kwargs.get('fixed_corners') is not None:
            Hi, ri = self.fixed_corners_constraints(**kwargs)
            H = sparse.vstack([H,Hi])
            r = np.hstack((r,ri))
            #print('2')
        if kwargs.get('self_closeness') != 0 and kwargs.get('self_closeness') is not None:
            Hi, ri = self.self_closeness_constraints(**kwargs)
            H = sparse.vstack([H,Hi])
            r = np.hstack((r,ri))
            #print('3')
        #print(kwargs.get('fixed_vertices'),kwargs.get('fixed_corners'),kwargs.get('self_closeness'))
        return H, r

    def iterative_constraints(self, **kwargs):
        N = kwargs.get('N')
        i = j = data = r = np.array([])
        H = sparse.coo_matrix((data,(i,j)), shape=(0,N))
        if kwargs.get('gliding'):
            Hi, ri = self.gliding_constraints(**kwargs)
            H = sparse.vstack((H, Hi))
            r = np.hstack((r, ri))
        if kwargs.get('reference_closeness'):
            Hi, ri = self.reference_closeness_constraints(**kwargs)
            H = sparse.vstack((H, Hi))
            r = np.hstack((r, ri)) 
        if kwargs.get('boundary_closeness'):
            Hi, ri = self.boundary_closeness_constraints(**kwargs)
            H = sparse.vstack((H, Hi))
            r = np.hstack((r, ri))   
        return H, r

    def fairness_energy(self, **kwargs):
        N = kwargs.get('N')
        i = j = data = s = np.array([])
        K = sparse.coo_matrix((data,(i,j)), shape=(0,N))
        if kwargs.get('mesh_fairness'):
            Ki, si = self.mesh_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
        if kwargs.get('tangential_fairness'):
            Ki, si = self.tangential_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
        if kwargs.get('boundary_fairness'):
            Ki, si = self.boundary_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
        if kwargs.get('spring_fairness'):
            Ki, si = self.spring_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
        if kwargs.get('corner_fairness'): # Hui add
            Ki, si = self.con_corner_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))  
        if kwargs.get('fairness_diagmesh'): # Hui add
            Ki, si = self.con_diagonal_mesh_fairness(**kwargs)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))               
        return K, s

    # -------------------------------------------------------------------------
    #                                 Relax
    # -------------------------------------------------------------------------

    def relax(self, iterations=1):
        self.glide(self.boundary_vertices())
        K, s = self.mesh_fairness(w=0.3)
        Ki, si = self.boundary_fairness(w=0.3)
        K = sparse.vstack((K, Ki))
        s = np.hstack((s, si))
        Ki, si = self.fixed_corners_constraints(w=2)
        K = sparse.vstack((K, Ki))
        s = np.hstack((s, si))
        for i in range(iterations):
            Ki, si = self.gliding_constraints(w=1)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
            Ki, si = self.reference_closeness_constraints(w=0.5)
            K = sparse.vstack((K, Ki))
            s = np.hstack((s, si))
            A = (K.T).dot(K)
            b = (K.T).dot(s)
            X = spsolve(A, b)
            V = self.V
            self.vertices[:,0] = X[0:V]
            self.vertices[:,1] = X[V:2*V]
            self.vertices[:,2] = X[2*V:3*V]

    def smooth_vertex_data(self, vertex_data, exclude_boundary=True, smoothing=1):
        "Hui note: all_data at each vertex-->[all_data; inner_data]"
        v_data = np.array(vertex_data)
        v, vj, l = self.vertex_ring_vertices_iterators(return_lengths=True)
        i = np.hstack((v, np.arange(self.V)))
        j = np.hstack((vj, np.arange(self.V)))
        data = np.hstack((np.ones(len(v)), -l)) * smoothing
        H1 = sparse.coo_matrix((data,(i,j)), shape=(self.V,self.V))
        vi = self.inner_vertices()
        W = len(vi)
        i = np.arange(W)
        data = np.ones(W)
        H2 = sparse.coo_matrix((data,(i,vi)), shape=(W, self.V))
        H = sparse.vstack((H1, H2))
        r = np.hstack((np.zeros(self.V), v_data[vi]))
        H = H.tocsc()
        A = (H.T).dot(H)
        b = (H.T).dot(r)
        X = spsolve(A, b)
        return X
    
    def smooth_quadface_data(self, face_data, rr=True,exclude_boundary=True, smoothing=1):
        """Hui add: similar to above smooth_vertex_data
        all faces, inner: rr_quadface\boundary_quadface
        """
        f_data = np.array(face_data)
        f, fj, l = self.face_ring_faces_iterators(return_lengths=True)
        
        i = np.hstack((f, np.arange(self.F)))
        j = np.hstack((fj, np.arange(self.F)))
        data = np.hstack((np.ones(len(f)), -l)) * smoothing
        H1 = sparse.coo_matrix((data,(i,j)), shape=(self.F,self.F))
        
        if rr:
            fi = self.rr_quadface_order
            
        if exclude_boundary:
            inn,_ = self.get_rr_quadface_boundaryquad_index()
            f_data = f_data[inn]
            fi = fi[inn]
            
        W = len(fi)
        i = np.arange(W)
        data = np.ones(W)
        H2 = sparse.coo_matrix((data,(i,fi)), shape=(W, self.F))
        
        H = sparse.vstack((H1, H2))
        r = np.hstack((np.zeros(self.F), f_data))

        H = H.tocsc()
        A = (H.T).dot(H)
        b = (H.T).dot(r)
        X = spsolve(A, b)
        return X

    def vertex_data_smoothing_matrix(self, offset=0, N=None, w=0.1):
        if N == None:
            N = self.V
        v, vj, l = self.vertex_ring_vertices_iterators(return_lengths=True)
        i = np.hstack((v, np.arange(self.V)))
        j = np.hstack((vj, np.arange(self.V))) + offset
        data = np.hstack((np.ones(len(v)), -l)) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(self.V, N))
        return H

    # -------------------------------------------------------------------------
    #                                 Remesh
    # -------------------------------------------------------------------------

    def incremental_remesh(self, target_length=1, iterations=10, relative=True):
        print('*** incremental remeshing ***')
        print(self)
        if relative:
            target_length = 0.8 * target_length * self.mean_edge_length()
        out = ('target length = {:.4f}\n').format(target_length)
        print(out)
        max_length = 4/3 * target_length
        min_length = 4/5 * target_length
        #reference_mesh = self.copy_mesh()
        self.set_reference()
        for i in range(iterations):
            self.split_edges(max_length)
            self.collapse_edges(min_length)
            self.equalize_valences()
            self.relax()
            print('iteration = ' + str(i+1))
        #self.reference_mesh = reference_mesh
        self.set_restart()
        print(self)
        print('*** done ***')

    # -------------------------------------------------------------------------
    #                                 Stress
    # -------------------------------------------------------------------------

    def principal_stresses(self):
        v, ej = self.vertex_ring_edges_iterators(sort=True)
        v, vj = self.vertex_ring_vertices_iterators(sort=True)
        Vi = self.vertices[v]
        Vj = self.vertices[vj]
        Vij = Vi - Vj
        W = self.force_densities
        W = W[ej]
        Wij = np.array([W, W, W]).T
        S = np.einsum('ij,ik -> ijk',Wij*Vij, Vij)
        S = utilities.sum_repeated(S,v)
        eig = np.linalg.eigh(S)
        Kmin = np.argmin(np.abs(eig[0]), axis=1)
        i = np.arange(self.V)
        i1 = np.array([2,0,1])[Kmin]
        i2 = np.array([1,2,0])[Kmin]
        S1 = eig[0][i,i1]
        S2 = eig[0][i,i2]
        V1x = eig[1][i,0,i1]
        V1y = eig[1][i,1,i1]
        V1z = eig[1][i,2,i1]
        V2x = eig[1][i,0,i2]
        V2y = eig[1][i,1,i2]
        V2z = eig[1][i,2,i2]
        V1 = np.array([V1x,V1y,V1z]).T
        V2 = np.array([V2x,V2y,V2z]).T
        N = self.vertex_normals()
        V1 = np.cross(N, V1)
        V1 = V1 / np.linalg.norm(V1, axis=1, keepdims=True)
        V2 = np.cross(N, V2)
        V2 = V2 / np.linalg.norm(V2, axis=1, keepdims=True)
        return (S1, S2, V2, V1)

    def make_crossfield_obj_file(self, V1, V2, file_name, scale=1):
        V = self.V
        mean = self.mean_edge_length() / 2
        V1 = np.reshape(V1,(V,3),order='F')
        V2 = np.reshape(V2,(V,3),order='F')
        P1a = self.vertices - scale * mean * V1
        P1b = self.vertices + scale * mean * V1
        P2a = self.vertices - scale * mean * V2
        P2b = self.vertices + scale * mean * V2
        P = np.vstack((P1a,P1b,P2a,P2b))
        name = ('{}.obj').format(file_name)
        obj = open(name, 'w')
        line = ('o {}\n').format('lines')
        obj.write((line))
        for v in range(4*V):
            vi = P[v]
            line = ('v {} {} {}\n').format(vi[0], vi[1], vi[2])
            obj.write((line))
        for v in range(V):
            line = ('l {} {}\n').format(v+1, v+V+1)
            obj.write((line))
        for v in range(V):
            line = ('l {} {}\n').format(2*V+v+1, v+3*V+1)
            obj.write((line))
        obj.close()

    # -------------------------------------------------------------------------
    #                               Minimal surface
    # -------------------------------------------------------------------------

    def minimize_area(self, threshold=1e-5):
        V = self.V
        v = self.boundary_vertices()
        W = len(v)
        j = np.hstack((v, V+v, 2*V+v))
        i = np.arange(W)
        i = np.hstack((i, W+i, 2*W+i))
        r = self.vertices[v,:]
        r = np.reshape(r, 3*W, order='F') * 10
        data = np.ones([3*W]) * 10
        H = sparse.coo_matrix((data, (i,j)), shape=(3*W, 3*V))
        K, s = self.mesh_fairness()
        A = sparse.vstack((K, H))
        b = np.hstack((s, r))
        A = A.tocsc()
        M = (A.T).dot(A)
        a = (A.T).dot(b)
        X  = spsolve(M, a)
        self.vertices[:,0] = X[0:V]
        self.vertices[:,1] = X[V:2*V]
        self.vertices[:,2] = X[2*V:3*V]
        self.mean_curvature_descent(threshold)

    def mean_curvature_descent(self, threshold=1e-5):
        A = self.area()
        M = self.mean_curvature_normal()
        eps = np.sum(np.linalg.norm(M, axis=1))
        i = 0
        lam = 1
        while True:
            stop = False
            while not stop:
                self.vertices += lam * M
                Mi = self.mean_curvature_normal()
                epsi = np.sum(np.linalg.norm(Mi, axis=1))
                if epsi < eps:
                    stop = True
                    M = Mi
                    #print((eps - epsi)/A)
                    if (eps - epsi)/A < threshold:
                        print(i)
                        return
                    eps = epsi
                else:
                    self.vertices -= lam * M
                    lam = 0.75*lam
                if lam < 1e-10:
                    return
            i += 1

    # -------------------------------------------------------------------------
    #                                 Dual
    # -------------------------------------------------------------------------

    def force_dual2(self):
        dual = self.dual_mesh()
        V = dual.V
        E = self.E
        H = self.halfedges
        mask = self.are_boundary_edges()
        I = E - np.sum(mask)
        mask0 = np.invert(mask)
        e, ind = np.unique(H[:,5], return_index=True)
        mask = mask0[H[ind,5]]
        ind = ind[mask]
        vectors = self.vertices[H[H[:,4],0]] - self.vertices[H[:,0]]
        versors = (vectors/np.linalg.norm(vectors,axis=1,keepdims=True))[ind]
        Hd = dual.halfedges
        i = np.arange(I)
        i = np.tile(i, 2)
        i = np.hstack((i, i+I, i+2*I))
        v1 = Hd[:,0]
        v2 = Hd[H[:,4],0]
        j = np.hstack((v1[ind], v2[ind]))
        j = np.hstack((j, V+j, 2*V+j))
        data = np.tile(np.hstack((np.ones(I), -np.ones(I))), 3)
        N = self.force_densities * self.edge_lengths()
        N = N[mask0]
        b = np.array([N]).T * versors
        b = np.reshape(b, (3*I), order='F')
        A = sparse.coo_matrix((data,(i,j)), shape=(3*I, 3*V))
        X = sparse.linalg.lsqr(A,b) ; X = X[0]
        dual.vertices = np.reshape(X, (V,3), order='F')
        return dual

    def force_dual(self, cut=True):
        primal = self.copy_gridshell()
        if cut:
            primal.cut()
        dual = primal.dual_mesh()
        V = dual.V
        E = primal.E
        H = primal.halfedges
        e, ind = np.unique(H[:,5], return_index=True)
        vectors = primal.vertices[H[H[:,4],0]] - primal.vertices[H[:,0]]
        vectors = vectors[ind]
        Hd = dual.halfedges
        i = np.arange(E)
        i = np.tile(i, 2)
        i = np.hstack((i, i+E, i+2*E))
        v1 = Hd[:,0]
        v2 = Hd[H[:,4],0]
        j = np.hstack((v1[ind], v2[ind]))
        j = np.hstack((j, V+j, 2*V+j))
        data = np.tile(np.hstack((np.ones(E), -np.ones(E))), 3)
        N = primal.force_densities
        b = np.array([N]).T * vectors
        b = np.reshape(b, (3*E), order='F')
        A = sparse.coo_matrix((data,(i,j)), shape=(3*E, 3*V))
        X = sparse.linalg.lsqr(A,b) ; X = X[0]
        dual.vertices = np.reshape(X, (V,3), order='F')
        return dual

    def force_dual3(self):
        primal = self.copy_gridshell()
        primal.cut()
        dual = primal.dual_mesh()
        V = dual.V
        E = primal.E
        H = primal.halfedges
        e, ind = np.unique(H[:,5], return_index=True)
        vectors = primal.vertices[H[H[:,4],0]] - primal.vertices[H[:,0]]
        versors = (vectors/np.linalg.norm(vectors,axis=1,keepdims=True))[ind]
        Hd = dual.halfedges
        i = np.arange(E)
        i = np.tile(i, 2)
        i = np.hstack((i, i+E, i+2*E))
        v1 = Hd[:,0]
        v2 = Hd[H[:,4],0]
        j = np.hstack((v1[ind], v2[ind]))
        j = np.hstack((j, V+j, 2*V+j))
        data = np.tile(np.hstack((np.ones(E), -np.ones(E))), 3)
        N = primal.force_densities * primal.edge_lengths()
        b = np.array([N]).T * versors
        b = np.reshape(b, (3*E), order='F')
        #t0 = time.time()
        A = sparse.coo_matrix((data,(i,j)), shape=(3*E, 3*V))
        '''
        A = A.tocsc()
        mask = self.are_boundary_edges()
        mask = np.invert(mask)
        mask = np.hstack((mask,mask,mask))
        A = A[mask]
        b = b[mask]
        M = (A.T).dot(A)
        a = (A.T).dot(b)
        X  = spsolve(M,a)
        '''
        X = sparse.linalg.lsqr(A,b) ; X = X[0]
        #print time.time()-t0
        dual.vertices = np.reshape(X, (V,3), order='F')
        return dual

    def cut2(self):
        H = self.halfedges
        densities = self.force_densities
        hc = super(GridshellNew, self).cut()
        if len(hc) == 0:
            return
        n_den = densities[H[hc,5]][:,0]
        self._force_densities = np.hstack((self._force_densities, n_den))

    # -------------------------------------------------------------------------
    #                                Connectivity
    # -------------------------------------------------------------------------

    def edge_side_edges(self):
        H = self.halfedges
        e, h = np.unique(H[:,5], return_index=True)
        e1 = H[H[H[h,2],2],5]
        e2 = H[H[H[H[h,4],2],2],5]
        e11 = H[H[H[H[H[H[h,2],2],4],2],2],5]
        e22 = H[H[H[H[H[H[H[h,4],2],2],4],2],2],5]
        A = H[h,1] != -1
        B = H[H[h,4],1] != -1
        mask = np.logical_and(A, B)
        e = e[mask]
        e1 = e1[mask]
        e2 = e2[mask]
        e11 = e11[mask]
        e22 = e22[mask]
        b = np.invert(self.are_boundary_edges())
        A = np.logical_and(b[e1], b[e2])
        B = np.logical_and(b[e11], b[e22])
        mask = np.logical_and(A, B)
        e = e[mask]
        e1 = e1[mask]
        e2 = e2[mask]
        e11 = e11[mask]
        e22 = e22[mask]
        return e, e1, e2, e11, e22

    def double_boundary_vertices(self):
        b = self.boundary_vertices()
        vi, vj = self.vertex_ring_vertices_iterators(sort=True)
        ring = np.copy(b)
        for v in b:
            ring = np.hstack((ring, (vj[vi==v])))
        return np.unique(ring)

