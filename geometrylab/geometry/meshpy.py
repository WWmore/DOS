#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import copy

from io import open

#import os

import numpy as np

from scipy.sparse import coo_matrix

from scipy import spatial

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.geometry.polyline import Polyline

from geometrylab.geometry.frame import Frame

from geometrylab.geometry.circle import circle_three_points

# -----------------------------------------------------------------------------

'''meshpy.py: The half-edge mesh data structure'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Half-Edge Mesh
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

class Mesh(object):

    def __init__(self, file_name=None):

        self._name = 'mesh'

        self._V = 0

        self._F = 0

        self._E = 0

        self._halfedges = None

        self._vertices = None

        self._uv = None

        self._corner_tolerance = 0.3

        self._kdtree = None

        self.__Htrash = None

        self.__Vtrash = None

        self.__Ftrash = None

        self.__Etrash = None

        if file_name is not None:
            self.read_obj_file(file_name)

    # -------------------------------------------------------------------------
    #                               Properties
    # -------------------------------------------------------------------------

    @property
    def type(self):
        return 'Mesh'

    @property
    def V(self):
        return self._V

    @property
    def F(self):
        return self._F

    @property
    def E(self):
        return self._E

    @property
    def vertices(self):
        return self._vertices

    @vertices.setter
    def vertices(self, vertices):
        vertices = np.array(vertices, 'f')
        if len(vertices.shape) != 2 or vertices.shape[1] != 3:
            out = '*vertices* attribute must be a np.ndarray of shape (n,3)!'
            raise ValueError(out)
        try:
            self._vertices[:,:] = vertices[:,:]
            self.geometry_update()
        except ValueError:
            self._vertices = vertices
            self.topology_update()

    @property
    def halfedges(self):
        return self._halfedges

    @halfedges.setter
    def halfedges(self, halfedges):
        halfedges = np.array(halfedges, 'i')
        if len(halfedges.shape) != 2 or halfedges.shape[1] != 6:
            out = '*halfedges* attribute must be a np.ndarray of shape (n,6)!'
            raise ValueError(out)
        self._halfedges = halfedges
        self.topology_update()

    @property
    def corner_tolerance(self):
        return self._corner_tolerance

    @corner_tolerance.setter
    def corner_tolerance(self, corner_tolerance):
        try:
            corner_tolerance = float(corner_tolerance)
        except:
            raise ValueError('*corner_tolerance* attribute must be a float!')
        if type(corner_tolerance) is float or type(corner_tolerance) is int:
            self._corner_tolerance = float(corner_tolerance)
        else:
            raise ValueError('*corner_tolerance* attribute must be a float!')

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if type(name) is str:
            self._name = name
        else:
            raise ValueError('*name* attribute must be a string!')

    @property
    def kdtree(self):
        if self._kdtree is None:
            self.make_kdtree()
        return self._kdtree

    @property
    def x(self):
        return self._vertices[:,0]


    @property
    def y(self):
        return self._vertices[:,1]

    @property
    def z(self):
        return self._vertices[:,2]

    # -------------------------------------------------------------------------
    #                               Updates
    # -------------------------------------------------------------------------

    def update_dimensions(self):
        self._V = np.amax(self.halfedges[:,0] + 1)
        self._F = np.amax(self.halfedges[:,1] + 1)
        self._E = np.amax(self.halfedges[:,5] + 1)

    def topology_update(self):
        self.update_dimensions()
        self._kdtree = None

    def geometry_update(self):
        self._kdtree = None

    # -------------------------------------------------------------------------
    #                               Methods
    # -------------------------------------------------------------------------

    def __str__(self):
        out = 'Halfedge mesh: |V| = {}, |F| = {}, |E| = {}'
        out = out.format(self.V, self.F, self.E)
        return out

    def copy_mesh(self):
        copy_mesh = Mesh()
        copy_mesh.name = self.name + '_copy'
        copy_mesh.corner_tolerance = self.corner_tolerance
        copy_mesh._V = self._V
        copy_mesh._F = self._F
        copy_mesh._E = self._E
        copy_mesh._vertices = np.copy(self.vertices)
        copy_mesh._halfedges = np.copy(self.halfedges)
        copy_mesh.__Htrash = copy.deepcopy(self.__Htrash)
        copy_mesh.__Vtrash = copy.deepcopy(self.__Vtrash)
        copy_mesh.__Ftrash = copy.deepcopy(self.__Ftrash)
        copy_mesh.__Etrash = copy.deepcopy(self.__Etrash)
        return copy_mesh

    # -------------------------------------------------------------------------
    #                          Build Data Structure
    # -------------------------------------------------------------------------

    def make_mesh(self, vertices_list, faces_list):
        try:
            self._make_halfedges(vertices_list, faces_list)
        except:
            try:
                faces_list = self.orient_faces(vertices_list, faces_list)
                self._make_halfedges(vertices_list, faces_list)
                print('*** Mesh reoriented ***')
            except:
                raise ValueError('The mesh is not manifold and orientable!')
        print(self)

    def _make_halfedges(self, vertices_list, faces_list):
        self._V = len(vertices_list)
        self._F = len(faces_list)
        self._vertices = np.array(vertices_list, 'f')
        orig = []
        face = []
        nexx = []
        prev = []
        twin_i = []
        twin_j = []
        h = 0
        for f in range(self.F):
            N = len(faces_list[f])
            orig.append(faces_list[f][0])
            face.append(f)
            nexx.append(h + 1)
            prev.append(h + N - 1)
            twin_i.append(faces_list[f][1])
            twin_j.append(faces_list[f][0])
            for v in range(1, N-1):
                orig.append(faces_list[f][v])
                face.append(f)
                nexx.append(h + v + 1)
                prev.append(h + v - 1)
                twin_i.append(faces_list[f][v+1])
                twin_j.append(faces_list[f][v])
            orig.append(faces_list[f][N-1])
            face.append(f)
            nexx.append(h)
            prev.append(h + N - 2)
            twin_i.append(faces_list[f][0])
            twin_j.append(faces_list[f][N-1])
            h += N
        H = np.zeros((h, 6), 'i')
        H[:,0] = orig
        H[:,1] = face
        H[:,2] = nexx
        H[:,3] = prev
        twin = coo_matrix((np.arange(h) + 1, (twin_i, twin_j)), shape=(h, h))
        twin = twin.tocsc()
        H[:,4] = twin[H[:,0],H[H[:,2],0]] - 1
        b = np.where(H[:,4] == -1)[0]
        boundary = H[b,:]
        boundary[:,0] = H[H[b,2],0]
        boundary[:,1] = -1
        boundary[:,4] = b
        #B = boundary.shape[0]  # test by Davide
        B = len(boundary)
        if B > 0:
            Bh = np.arange(h, h+B)
            H[b,4] = Bh
            zeros = np.zeros(B)
            p = coo_matrix((Bh, (H[b,0], zeros)), shape=(self.V, 1))
            p = p.tocsc()
            # print(boundary[:,0]) # test by Davide
            # print(zeros) # test by Davide
            # print(p) # test by Davide
            boundary[:,3] = p[boundary[:,0],zeros]
            i = boundary[boundary[:,3]-h,0]
            p = coo_matrix((Bh, (i, zeros)), shape=(self.V, 1))
            p = p.tocsc()
            boundary[:,2] = p[boundary[:,0], zeros]
            H = np.vstack((H, boundary))
        K = H[:,(3,4)]
        K[:,0] = np.arange(H.shape[0])
        m = np.amin(K, axis=1)
        u = np.unique(m)
        imap = np.arange(np.max(u) + 1, dtype=int)
        imap[u] = np.arange(u.shape[0])
        H[:,5] = imap[m]
        self._halfedges = H
        self._E = int(H.shape[0] / 2)
        self.topology_update()

    # -------------------------------------------------------------------------
    #                                Navigating
    # -------------------------------------------------------------------------

    def origin(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,0]
        return H[halfedge_index,0]

    def face(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,1]
        return H[halfedge_index,1]

    def next(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,2]
        return H[halfedge_index,2]

    def previous(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,3]
        return H[halfedge_index,3]

    def twin(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,4]
        return H[halfedge_index,4]

    def edge(self, halfedge_index=None):
        H = self.halfedges
        if halfedge_index is None:
            return H[:,5]
        return H[halfedge_index,5]

    # -------------------------------------------------------------------------
    #                               Reading
    # -------------------------------------------------------------------------

    def read_obj_file(self, file_name):
        file_name = str(file_name)
        self.name = file_name.split('.')[0]
        obj_file = open(file_name, encoding='utf-8')
        vertices_list = []
        uv_list = []
        faces_list = []
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
            elif splited_line[0] == 'f':
                v_list = []
                L = len(splited_line)
                try:
                    for i in range(1, L):
                        splited_face_data = splited_line[i].split('/')
                        v_list.append(int(splited_face_data[0]) - 1 )
                    faces_list.append(v_list)
                except ValueError:
                    v_list = []
                    for i in range(1, L-1):
                        v_list.append(int(splited_line[i]) - 1 )
                    faces_list.append(v_list)
            if splited_line[0] == 'vt':
                split_u = splited_line[1].split('\n')
                u = float(split_u[0])
                split_v = splited_line[2].split('\n')
                v = float(split_v[0])
                vertices_list.append([u,v])
            if len(uv_list) > 0:
                self._uv = np.array(uv_list)
        self.make_mesh(vertices_list, faces_list)

    # -------------------------------------------------------------------------
    #                               Writing
    # -------------------------------------------------------------------------

    def make_obj_file(self, file_name, overwrite=False):
        path = utilities.make_filepath(file_name, 'obj', overwrite)
        obj = open(path, 'w')
        line = ('o {}\n').format(file_name)
        obj.write(line)
        faces = self.faces_list()
        for v in range(self.V):
            vi = self.vertices[v]
            line = ('v {} {} {}\n').format(vi[0], vi[1], vi[2])
            obj.write(line)
        for f in range(self.F):
             obj.write('f ')
             N = len(faces[f])
             for v in range(N - 1):
                 vf = str(faces[f][v] + 1)
                 #obj.write(unicode(vf + '//' + ' '))
                 obj.write(vf + ' ')
             vf = str(faces[f][N - 1] + 1)
             #obj.write(unicode(vf + '//' + '\n'))
             obj.write(vf + '\n')
        obj.close()
        return path.split('.')[0]

    # -------------------------------------------------------------------------
    #                                 Ordering
    # -------------------------------------------------------------------------

    def vertex_ring_ordered_halfedges(self):
        H = np.copy(self.halfedges)
        i = np.argsort(H[:,0])
        v = H[i,0]
        index = np.arange(H.shape[0])
        _, j = np.unique(v, True)
        v = np.delete(v,j)
        index = np.delete(index,j)
        while v.shape[0] > 0:
            _, j = np.unique(v, True)
            i[index[j]] = H[H[i[index[j] - 1],3],4]
            v = np.delete(v,j)
            index = np.delete(index,j)
        return i

    def face_ordered_halfedges(self):
        H = np.copy(self.halfedges)
        i = np.argsort(H[:,1])
        i = i[np.where(H[i,1] >= 0)]
        f = H[i,1]
        index = np.arange(i.shape[0])
        _, j = np.unique(f, True)
        f = np.delete(f,j)
        index = np.delete(index, j)
        while f.shape[0] > 0:
            _, j = np.unique(f, True)
            i[index[j]] = H[i[index[j] - 1],2]
            f = np.delete(f, j)
            index = np.delete(index, j)
        return i

    # -------------------------------------------------------------------------
    #                              Iterators
    # -------------------------------------------------------------------------

    def vertex_ring_vertices_iterators(self, sort=False, order=False,
                                       return_lengths=False):
        H = self.halfedges
        v  = H[:,0]
        vj = H[H[:,4],0]
        if order:
            i  = self.vertex_ring_ordered_halfedges()
            v  = v[i]
            vj = vj[i]
        elif sort:
            i  = np.argsort(v)
            v  = v[i]
            vj = vj[i]
        if return_lengths:
            i  = np.ones(vj.shape[0], dtype=int)
            lj = utilities.sum_repeated(i,v)
            return v, vj, lj
        else:
            return v, vj

    def vertex_ring_faces_iterators(self, sort=False, order=False):
        H = self.halfedges
        if order:
            i  = self.vertex_ring_ordered_halfedges()
            v  = H[i,0]
            fj = H[i,1]
        else:
            i  = np.where(H[:,1] >= 0)[0]
            v  = H[i,0]
            fj = H[i,1]
            if sort:
                i  = np.argsort(v)
                v  = v[i]
                fj = fj[i]
        return v, fj

    def vertex_ring_edges_iterators(self, sort=False, order=False):
        H = self.halfedges
        v  = H[:,0]
        ej = H[:,5]
        if order:
            i  = self.vertex_ring_ordered_halfedges()
            v  = v[i]
            ej = ej[i]
        elif sort:
            i  = np.argsort(v)
            v  = v[i]
            ej = ej[i]
        return v, ej

    def face_edge_vertices_iterators(self, sort=False, order=False):
        H = self.halfedges
        f  = H[:,1]
        vi = H[:,0]
        vj = H[H[:,2],0]
        if order:
            i  = self.face_ordered_halfedges()
            f  = f[i]
            vi = vi[i]
            vj = vj[i]
        else:
            i  = np.where(H[:,1] >= 0)[0]
            f  = f[i]
            vi = vi[i]
            vj = vj[i]
            if sort:
               i  = np.argsort(f)
               vi = vi[i]
               vj = vj[i]
        return f, vi, vj

    def face_vertices_iterators(self):
        H = self.halfedges
        i  = self.face_ordered_halfedges()
        vi = H[i,0]
        f  = H[i,1]
        return f, vi

    def face_edges_iterators(self):
        H = self.halfedges
        i  = self.face_ordered_halfedges()
        ei = H[i,5]
        f  = H[i,1]
        return f, ei

    def edge_vertices_iterators(self, sort=False):
        H = self.halfedges
        e  = H[:,5]
        vi = H[:,0]
        if sort:
            i  = np.argsort(H[:,5])
            e  = e[i]
            vi = vi[i]
        return e, vi

    def vertex_double_ring_vertices_iterators(self):
        #import time
        #t0 = time.time()
        v, vj = self.vertex_ring_vertices_iterators(sort=True)
        M = coo_matrix((vj, (v, vj)), shape=(self.V, self.V))
        M = M.todense()
        ring = np.copy(M)
        while v.shape[0] > 0:
            vi, j = np.unique(v, True)
            ring[vi] += M[vj[j]]
            v = np.delete(v, j)
            vj = np.delete(vj, j)
        #t4 = time.time()
        #print(t4-t0)
        return ring.nonzero()

    # -------------------------------------------------------------------------
    #                               Ring lists
    # -------------------------------------------------------------------------

    def vertex_ring_vertices_list(self):
        ring_list = [[] for i in range(self.V)]
        v, vj = self.vertex_ring_vertices_iterators(order=True)
        for i in range(len(v)):
            ring_list[v[i]].append(vj[i])
        return ring_list

    def vertex_double_ring_vertices_list(self):
        ring_list = [[] for i in range(self.V)]
        v, vj = self.vertex_double_ring_vertices_iterators()
        for i in range(len(v)):
            ring_list[v[i]].append(vj[i])
        return ring_list

    def vertex_ring_edges_list(self):
        ring_list = [[] for i in range(self.V)]
        v, ej = self.vertex_ring_edges_iterators(order=True)
        for i in range(len(v)):
            ring_list[v[i]].append(ej[i])
        return ring_list

    def vertex_ring_faces_list(self):
        ring_list = [[] for i in range(self.V)]
        v, fj = self.vertex_ring_faces_iterators(order=True)
        for i in range(len(v)):
            ring_list[v[i]].append(fj[i])
        return ring_list

    #--------------------------------------------------------------------------
    #                                 Faces
    #--------------------------------------------------------------------------

    def face_lengths(self):
        H = self.halfedges
        f = H[H[:,1] >= 0,1]
        f = f[np.argsort(f)]
        i = np.ones((f.shape), 'i')
        lengths = utilities.sum_repeated(i, f)
        return lengths

    def cell_arrays(self):
        H = self.halfedges
        i  = self.face_ordered_halfedges()
        vi = H[i,0]
        f  = H[i,1]
        i = np.ones((f.shape[0]), 'i')
        j = np.arange(f.shape[0])
        _, k = np.unique(f, True)
        lengths = utilities.sum_repeated(i, f)
        index = j[k]
        cells = np.insert(vi, index, lengths)
        cell_types = lengths - 3
        cell_types[np.where(cell_types[:] > 2)[0]] = 2
        return cells, cell_types

    def faces_list(self):
        faces_list = [[] for i in range(self.F)]
        fi, vj = self.face_vertices_iterators()
        for i in range(len(fi)):
            faces_list[fi[i]].append(vj[i])
        return faces_list

    def face_triangles(self):
        H = np.copy(self.halfedges)
        h = np.argsort(H[:,1])
        h = h[np.where(H[h,1] >= 0)]
        f = H[h,1]
        f_i, j = np.unique(f, True)
        one = np.arange(j.shape[0])
        f = np.delete(f, j)
        f = np.delete(f, j-one)
        f = np.delete(f, j-2*one)
        T = np.column_stack((H[j,0], H[H[j,2],0], H[H[H[j,2],2],0]))
        nex = H[H[H[j,2],2],2]
        face_index = f_i
        offset = 0
        while len(f) > 0:
            f_i, j = np.unique(f, True)
            T_i = np.column_stack((T[offset+f_i,-1], H[nex[f_i],0], T[f_i,0]))
            f = np.delete(f, j)
            nex = H[nex,2]
            T = np.vstack((T, T_i))
            face_index = np.hstack((face_index, f_i))
            offset += len(f_i)
        return T, face_index

    # -------------------------------------------------------------------------
    #                                  Edges
    # -------------------------------------------------------------------------

    def edge_vertices(self):
        H  = self.halfedges
        v  = H[np.argsort(H[:,5]),0]
        v1 = v[0::2]
        v2 = v[1::2]
        return v1, v2

    def edge_faces(self):
        H  = self.halfedges
        f  = H[np.argsort(H[:,5]),1]
        f1 = f[0::2]
        f2 = f[1::2]
        return f1, f2

    def vertices_edge_map(self):
        H  = self.halfedges
        v1 = H[:,0]
        v2 = H[H[:,4],0]
        e  = H[:,5]
        edge_map = coo_matrix((e, (v1,v2)), shape=(self.V, self.V))
        edge_map = edge_map.tocsc()
        return edge_map

    def vertices_edge_faces_maps(self):
        H  = self.halfedges
        v1 = H[:,0]
        v2 = H[H[:,4],0]
        f1 = H[:,1]
        f2 = H[H[:,4],1]
        f1Map = coo_matrix((f1, (v1,v2)), shape=(self.V, self.V))
        f2Map = coo_matrix((f2, (v1,v2)), shape=(self.V, self.V))
        f1Map = f1Map.tocsc()
        f2Map = f2Map.tocsc()
        return f1Map, f2Map

    def edge_halfedges(self):
        H = self.halfedges
        e = np.argsort(H[:,5])
        h1 = e[0::2]
        h2 = e[1::2]
        return h1, h2

    # -------------------------------------------------------------------------
    #                                 Boundary
    # -------------------------------------------------------------------------

    def boundary_vertices(self):
        H = self.halfedges
        b = np.where(H[:,1] == -1)[0]
        v = H[b,0]
        return v

    def inner_halfedges(self):
        H = self.halfedges
        h = np.where(H[:,1] != -1)[0]
        return h

    def boundary_halfedges(self):
        H = self.halfedges
        b = np.where(H[:,1] == -1)[0]
        return b

    def boundary_faces(self):
        H = self.halfedges
        b = self.boundary_halfedges()
        e = H[H[b,4],1]
        return e

    def inner_vertices(self):
        b = self.boundary_vertices()
        v = np.arange(self.V)
        mask = np.invert(np.in1d(v, b))
        v = v[mask]
        return v

    def boundary_curves(self, corner_split=False):
        H = self.halfedges
        boundaries = []
        boundary_halfedges = []
        for h in range(H.shape[0]):
            if H[h,1] == -1 and h not in boundary_halfedges:
                boundary = []
                h_he = h
                boundary_halfedges.append(h_he)
                boundary.append(H[h_he,0])
                h_he = H[h_he,2]
                while h_he != h:
                    boundary_halfedges.append(h_he)
                    boundary.append(H[h_he,0])
                    h_he = H[h_he,2]
                boundaries.append(np.array(boundary))
        if corner_split:
            corner_boundaries = []
            corners = self.mesh_corners()
            for boundary in boundaries:
                indices = np.arange(len(boundary))
                c = indices[np.in1d(boundary, corners)]
                boundary = np.split(boundary, c)
                for i in range(len(boundary) - 1):
                    a = boundary[i]
                    boundary[i] = np.insert(a, a.shape ,boundary[i+1][0])
                if len(boundary) > 1:
                    boundary[0] = np.hstack((boundary[-1], boundary[0]))
                    del boundary[-1]
                corner_boundaries.extend(boundary)
            boundaries = corner_boundaries
        return boundaries

    def boundary_curves_halfedges(self, corner_split=False):
        H = self.halfedges
        boundaries = []
        visited = []
        for h in range(H.shape[0]):
            if H[h,1] == -1 and h not in visited:
                boundary = []
                h_he = h
                boundary.append(h_he)
                h_he = H[h_he,2]
                while h_he != h:
                    boundary.append(h_he)
                    h_he = H[h_he,2]
                boundaries.append(np.array(boundary))
                visited.extend(boundary)
        if corner_split:
            corner_boundaries = []
            corners = self.mesh_corners()
            for boundary in boundaries:
                indices = np.arange(len(boundary))
                c = indices[np.in1d(H[boundary,0], corners)]
                boundary = np.split(boundary, c)
                if len(boundary) > 1:
                    boundary[0] = np.hstack((boundary[-1], boundary[0]))
                    del boundary[-1]
                corner_boundaries.extend(boundary)
            boundaries = corner_boundaries
        return boundaries

    def boundary_polylines(self):
        polylines = []
        curves = self.boundary_curves(corner_split=False)
        for curve in curves:
            polyline = Polyline(self.vertices[curve,:], closed=True)
            polyline.corner_tolerance = self.corner_tolerance
            polylines.append(polyline)
        return polylines

    def mesh_corners(self):
        H = self.halfedges
        b = np.where(H[:,1] == -1)[0]
        v0 = H[b,0]
        vp = H[H[b,3],0]
        vn = H[H[b,2],0]
        Vp = self.vertices[v0,:] - self.vertices[vp,:]
        Vn = self.vertices[vn,:] - self.vertices[v0,:]
        Vp = Vp / np.linalg.norm(Vp, axis=1, keepdims=True)
        Vn = Vn / np.linalg.norm(Vn, axis=1, keepdims=True)
        C = np.einsum('ij,ij->i', Vp, Vn)
        corners = v0[np.where(C[:] < self.corner_tolerance)[0]]
        return corners

    def double_boundary_vertices(self):
        bf = self.are_boundary_faces()
        v, fj = self.vertex_ring_faces_iterators(sort=True)
        bf = bf[fj]
        bf = utilities.sum_repeated(bf, v)
        bf = np.where(bf > 0)[0]
        return bf

    def boundary_edges(self):
        H = self.halfedges
        ind = np.where(H[:,1] == -1)
        e = H[ind,5]
        e = np.unique(e)
        return e

    # -------------------------------------------------------------------------
    #                              Global queries
    # -------------------------------------------------------------------------

    def are_boundary_edges(self):
        H = self.halfedges
        B = H[np.argsort(H[:,5])]
        B = B[:,1] == -1
        bound = np.logical_or(B[0::2], B[1::2])
        return bound

    def are_boundary_faces(self):
        H = self.halfedges
        f = np.where(H[:,1] != -1)[0]
        B = H[H[f,4],1] == -1
        i = np.argsort(H[f,1])
        bound = utilities.sum_repeated(B[i], H[f,1])
        return bound

    # -------------------------------------------------------------------------
    #                                Normals
    # -------------------------------------------------------------------------

    def face_vector_areas(self):
        f, v1, v2 = self.face_edge_vertices_iterators(order=True)
        V1 = self.vertices[v1,:]
        V2 = self.vertices[v2,:]
        N  = np.cross(V1,V2)
        normals = utilities.sum_repeated(N, f)
        return 0.5 * normals

    def face_normals(self):
        N = self.face_vector_areas()
        N = N / np.linalg.norm(N, axis=1, keepdims=True)
        return N

    def vertex_normals(self):
        N = self.face_vector_areas()
        v, fi = self.vertex_ring_faces_iterators(sort=True)
        N = N[fi,:]
        normals = utilities.sum_repeated(N, v)
        normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)
        return normals

    def edge_normals(self):
        N = self.face_normals()
        N = np.insert(N, N.shape[0], 0, axis=0)
        f1, f2 = self.edge_faces()
        normals = N[f1] + N[f2]
        normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)
        return normals

    def boundary_normals(self):
        H = self.halfedges
        b = np.where(H[:,1] == -1)[0]
        face_normals = self.face_normals()
        N1 = face_normals[H[H[b,4],1]]
        N2 = face_normals[H[H[H[b,3],4],1]]
        normals = np.zeros((self.V, 3))
        E1 = self.vertices[H[H[b,2],0]] - self.vertices[H[b,0]]
        E2 = self.vertices[H[b,0]] - self.vertices[H[H[b,3],0]]
        N = np.cross(N1,E1) + np.cross(N2,E2)
        N = N / np.linalg.norm(N, axis=1, keepdims=True)
        normals[H[b,0],:] = N
        return normals

    def boundary_tangents(self, normalize=True):
        H = self.halfedges
        b = np.where(H[:,1] == -1)[0]
        V1 = self.vertices[H[H[b,3],0]]
        V2 = self.vertices[H[H[b,2],0]]
        T = (V2 - V1)
        if normalize:
            T = T / np.linalg.norm(T, keepdims=True)
        else:
            T = T/2
        tangents = np.zeros((self.V, 3))
        tangents[H[b,0],:] = T
        return tangents

    # -------------------------------------------------------------------------
    #                                  Area
    # -------------------------------------------------------------------------

    def face_areas(self):
        N = self.face_vector_areas()
        A = np.linalg.norm(N, axis=1)
        return A

    def area(self):
        A = self.face_areas()
        A = np.sum(A)
        return A

    def vertex_ring_areas(self):
        L = self.face_lengths()
        A = self.face_areas()
        v, fi = self.vertex_ring_faces_iterators(sort=True)
        ring_area = A[fi]/L[fi]
        ring_area = utilities.sum_repeated(ring_area, v)
        return ring_area

    # -------------------------------------------------------------------------
    #                               Closeness
    # -------------------------------------------------------------------------

    def make_kdtree(self):
        kdtree = spatial.cKDTree(self.vertices)
        self._kdtree = kdtree

    def closest_vertices(self, points, make_tree=False):
        if self._kdtree is None:
            self.make_kdtree()
        elif make_tree:
            self.make_kdtree()
        closest = self._kdtree.query(points)[1]
        return closest

    # -------------------------------------------------------------------------
    #                                Geometry
    # -------------------------------------------------------------------------

    def edge_mid_points(self):
        v1, v2 = self.edge_vertices()
        M = 0.5*(self.vertices[v1] + self.vertices[v2])
        return M

    def edge_vectors(self):
        v1, v2 = self.edge_vertices()
        Ev = self.vertices[v1] - self.vertices[v2]
        return Ev

    def face_barycenters(self):
        H = self.halfedges
        H = H[np.where(H[:,1] >= 0)[0],:]
        i = np.argsort(H[:,1])
        f = H[i,1]
        v = H[i,0]
        V = self.vertices[v,:]
        B = utilities.sum_repeated(V,f)
        L = self.face_lengths()
        L = np.column_stack((L,L,L))
        B = B/L
        return B

    def edge_lengths(self):
        v1, v2 = self.edge_vertices()
        V1 = self.vertices[v1]
        V2 = self.vertices[v2]
        V = V1 - V2
        L = np.linalg.norm(V, axis=1)
        return L

    def mean_edge_length(self):
        L = self.edge_lengths()
        M = np.sum(L) / L.shape[0]
        return M

    def face_planarity(self, scale_invariant=True):
        planarity = np.zeros((self.F))
        f, vi = self.face_vertices_iterators()
        i = np.ones((f.shape[0]),dtype=int)
        j = np.arange(f.shape[0])
        _, k = np.unique(f, True)
        L = utilities.sum_repeated(i, f)
        index = j[k]
        quad = np.where(L > 3)[0]
        shift = 0
        while len(quad) > 0:
            P1 = self.vertices[vi[index[quad] + shift]]
            P2 = self.vertices[vi[index[quad] + shift + 1]]
            P3 = self.vertices[vi[index[quad] + shift + 2]]
            P4 = self.vertices[vi[index[quad] + shift + 3]]
            V1 = P3 - P1
            V2 = P4 - P2
            N  = np.cross(V1,V2)
            eps = np.finfo(float).eps
            P12 = P2 - P1
            norm = ((np.linalg.norm(N, axis=1) + eps))
            d = np.einsum('ij,ij->i', P12, N) / norm
            if scale_invariant:
                d1 = np.linalg.norm(V1, axis=1)
                d2 = np.linalg.norm(V1, axis=1)
                p = np.abs(d) / ((d1 + d2)/2)
            else:
                p = np.abs(d)
            planarity[quad] = np.maximum(p, planarity[quad])
            L -= 1
            shift += 1
            quad = np.where(L > 3)[0]
        return planarity

    def edge_versors(self):
        v1, v2 = self.edge_vertices()
        ver = self.vertices[v2] - self.vertices[v1]
        ver = ver / np.linalg.norm(ver, axis=1, keepdims=True)
        return ver

    def bounding_box(self):
        Xmin = np.min(self.vertices[:,0])
        Xmax = np.max(self.vertices[:,0])
        Ymin = np.min(self.vertices[:,1])
        Ymax = np.max(self.vertices[:,1])
        Zmin = np.min(self.vertices[:,2])
        Zmax = np.max(self.vertices[:,2])
        return ([Xmin, Xmax], [Ymin, Ymax], [Zmin, Zmax])

    def mesh_center(self):
        B = self.bounding_box()
        x = (B[0][0] + B[0][1])/2
        y = (B[1][0] + B[1][1])/2
        z = (B[2][0] + B[2][1])/2
        return np.array([x,y,z])

    def face_circum_circles(self):
        f, vi = self.face_vertices_iterators()
        _, j = np.unique(f,True)
        p1 = self.vertices[vi[j]]
        p2 = self.vertices[vi[j+1]]
        p3 = self.vertices[vi[j+2]]
        circles = circle_three_points(p1, p2, p3)
        return circles

    # -------------------------------------------------------------------------
    #                                 Topology
    # -------------------------------------------------------------------------

    def flip_normals(self):
        H = self.halfedges
        H[:,0] = H[H[:,2],0]
        H[:,[2,3]] = H[:,[3,2]]

    def orient_faces(self, vertices_list, faces_list):
        F = len(faces_list)
        V = len(vertices_list)
        fmap = -np.ones((V,V), dtype=int)
        inconsistent = np.zeros((V,V), dtype=int)
        flipped = np.zeros(F, dtype=bool)
        oriented = np.zeros(F, dtype=bool)
        oriented_faces = copy.deepcopy(faces_list)
        for f in range(F):
            face = faces_list[f]
            for j in range(len(face)):
                v0 = face[j-1]
                v1 = face[j]
                if fmap[v0,v1] == -1:
                    fmap[v0,v1] = f
                else:
                    fmap[v1,v0] = f
                    inconsistent[v0,v1] = True
                    inconsistent[v1,v0] = True
        ring = [0]
        oriented[0] = True
        i = 1
        while len(ring) > 0:
            next_ring = []
            for f in ring:
                face = faces_list[f]
                for j in range(len(face)):
                    flip = False
                    v0 = face[j-1]
                    v1 = face[j]
                    if fmap[v0,v1] == f:
                        v2 = v1
                        v3 = v0
                    else:
                        v2 = v0
                        v3 = v1
                    if inconsistent[v2,v3] and not flipped[f]:
                        flip = True
                    if not inconsistent[v2,v3] and flipped[f]:
                        flip = True
                    fi = fmap[v2,v3]
                    if fi != -1 and not oriented[fi]:
                        if fi not in next_ring:
                            next_ring.append(fi)
                        if flip:
                            oriented_faces[fi].reverse()
                            flipped[fi] = True
                        i += 1
                        oriented[fi] = True
                        if i == F:
                            return oriented_faces
            ring = next_ring
            if len(ring) == 0:
                try:
                    ring = [np.where(oriented == False)[0][0]]
                except:
                    return

    def vertex_ring_expansion(self, v_index, callback=None, depth=None):
        vi, vj = self.vertex_ring_vertices_iterators()
        mring = np.full(self.V, False)
        sring = np.full(self.V, False)
        search = np.array([v_index], dtype='i')
        mring[v_index] = True
        if depth is None:
            depth = self.V
        for i in range(depth):
            sring[:] = False
            for v in search:
                ring = vj[vi == v]
                ring = ring[np.invert(mring[ring])]
                mring[ring] = True
                sring[ring] = True
                if callable(callback):
                    callback(v, ring)
            search = np.where(sring)[0]
            if np.all(mring):
                return np.where(mring)[0]
        return np.where(mring)[0]

    def make_simply_connected(self):
        curves = self.boundary_curves(corner_split=False)
        while len(curves) > 1:
            curve = curves[0]
            i = 0
            v = curve[i]
            while len(self.vertex_ring_vertices(v)) < 3:
                i += 1
                v = curve[i]
            self.cut(v)
            curves = self.boundary_curves(corner_split=False)

    # -------------------------------------------------------------------------
    #                                 Curves
    # -------------------------------------------------------------------------

    def mesh_curves(self):
        _,_, valence = self.vertex_ring_vertices_iterators(return_lengths=True)
        boundary_vertices = self.boundary_vertices()
        H = self.halfedges
        boundaries = self.boundary_curves_halfedges(True)
        done = []
        curves = []
        for boundary in boundaries:
            family = []
            for h in boundary:
                if H[h,0] not in done:
                    curve = [H[h,0]]
                    if valence[H[h,0]] <= 3:
                        turn = 1
                    else:
                        turn = 2
                    for i in range(turn):
                        h = H[H[h,4],2]
                    vertex = H[H[h,4],0]
                    stop = False
                    exclude = False
                    if vertex in boundary_vertices:
                        stop = True
                        exclude = True
                    while not stop:
                        curve.append(vertex)
                        if vertex in boundary_vertices:
                            stop = True
                            done.append(vertex)
                        if valence[vertex] <= 4:
                            turn = 1
                        else:
                            turn = 2
                        for i in range(turn):
                            h = H[H[H[h,2],4],2]
                        vertex = H[H[h,4],0]
                    if not exclude:
                        family.append(curve)
            if len(family) > 0:
                curves.append(family)
        curves.append(self.boundary_curves(True))
        return curves

    def mesh_polylines(self):
        curves = self.mesh_curves()
        polylines = []
        for family in curves:
            poly_family = []
            for curve in family:
                poly_family.append(Polyline(self.vertices[curve,:]))
            polylines.append(poly_family)
        return polylines

    # -------------------------------------------------------------------------
    #                             Transformations
    # -------------------------------------------------------------------------

    def move(self, displacement_vector):
        self.vertices[:,[0,1,2]] += np.array(displacement_vector)[[0,1,2]]

    def scale(self, factor, center=[0,0,0]):
        self.vertices[:,:] *= factor
        self.vertices[:,0] -= center[0]
        self.vertices[:,1] -= center[1]
        self.vertices[:,2] -= center[2]

    # -------------------------------------------------------------------------
    #                     Discrete Differential Geometry
    # -------------------------------------------------------------------------

    def vertex_ring_parametrization(self):
        v, vj, l = self.vertex_ring_vertices_iterators(sort=True,
                                                    return_lengths=True)
        index = np.arange(v.shape[0])
        step = np.zeros(v.shape[0])
        _, unique = np.unique(v, return_index=True)
        vertices = np.delete(v, unique)
        index = np.delete(index, unique)
        value = 0
        while len(vertices) > 0:
            value += 1
            _, unique = np.unique(vertices, return_index=True)
            step[index[unique]] = value
            vertices = np.delete(vertices, unique)
            index = np.delete(index, unique)
        phi = 2*np.pi*step / l[v]
        U = np.sin(phi)
        V = np.cos(phi)
        return v, vj, U, V

    def vertex_local_frame(self):
        o = self.vertices
        z = self.vertex_normals()
        x = utilities.orthogonal_vectors(z)
        y = np.cross(z,x)
        frame = Frame(o, x, y, z)
        return frame

    def edge_angle_vectors(self):
        sin = self.edge_sine_vectors()
        beta = np.arcsin(np.linalg.norm(sin, axis=1))
        beta = np.array([beta,beta,beta]).T * utilities.normalize(sin)
        return beta

    def edge_angles(self):
        sin = self.edge_sine_vectors()
        beta = np.arcsin(np.linalg.norm(sin, axis=1))
        return beta

    def edge_sine_vectors(self):
        v, ej = self.vertex_ring_edges_iterators(sort=True)
        normals = self.face_normals()
        v1, v2 = self.edge_vertices()
        f1, f2 = self.edge_faces()
        bf1 = np.where(f1 == -1)[0]
        bf2 = np.where(f2 == -1)[0]
        f1[bf1] = f2[bf1]
        f2[bf2] = f1[bf2]
        F1 = normals[f1,:]
        F2 = normals[f2,:]
        sin = np.cross(F1, F2)
        return sin

    def extended_shape_operator(self, area_normalization=False, use_sine=False):
        v, ej = self.vertex_ring_edges_iterators(sort=True)
        v1, v2 = self.edge_vertices()
        if use_sine:
            B = self.edge_sine_vectors()
        else:
            B = self.edge_angle_vectors()
        V12 = self.vertices[v1[ej]] - self.vertices[v2[ej]]
        B12 = B[ej] / 2
        W = np.einsum('ij,ik -> ijk', B12, V12)
        W = utilities.sum_repeated(W, v)
        if area_normalization:
            A = self.vertex_ring_areas()
            A = np.array([[A,A,A],[A,A,A],[A,A,A]]).T
            W = W/A
        return W

    def principal_curvatures(self, area_normalization=False, use_sine=False):
        W = self.extended_shape_operator(area_normalization, use_sine)
        try:
            eig = np.linalg.eigh(W)
            srt = np.argsort(np.abs(eig[0]), axis=1)
            i = np.arange(self.V)
            i1 = srt[:,1]
            i2 = srt[:,2]
            k1 = eig[0][i,i1]
            k2 = eig[0][i,i2]
            V1 = eig[1][i,:,i1]
            V2 = eig[1][i,:,i2]
            N = self.vertex_normals()
            D1 = utilities.normalize(np.cross(V1,N))
            D2 = utilities.normalize(np.cross(V2,N))
        except:
            V = self.V
            return (np.ones(V), np.ones(V), np.ones((V,3)), np.ones((V,3)))
        return (k1, k2, D1, D2)

    def curvature_ratios(self):
        K = self.principal_curvatures(True)
        k1 = K[0]
        k2 = K[1]
        R = np.maximum(np.abs(k2/(k1 + 1e-10)),
                               np.abs(k1/(k2 + 1e-10))) * np.sign(k1*k2)
        return R

    def gaussian_curvature(self):
        K = self.principal_curvatures(True, use_sine=True)
        return K[0]*K[1]

    def mean_curvature(self):
        K = self.principal_curvatures(True)
        H =  0.5*(K[0] + K[1])
        return H

    def edge_cotangents_weigths(self):
        H = self.halfedges
        V0 = self.vertices[H[:,0]]
        V1 = self.vertices[H[H[:,2],0]]
        V2 = self.vertices[H[H[:,3],0]]
        E1 = (V0 - V2) / np.linalg.norm(V0 - V2, axis=1, keepdims=True)
        E2 = (V1 - V2) / np.linalg.norm(V1 - V2, axis=1, keepdims=True)
        cos = np.einsum('ij,ij->i', E1, E2)
        sin = np.linalg.norm(np.cross(E1, E2, axis=1), axis=1)
        cotan = cos / sin
        h = np.argsort(H[:,5])
        cotan = cotan[h]
        e = H[h,5]
        cotan = utilities.sum_repeated(cotan, e)
        return cotan / 2

    def mean_curvature_normal(self):
        b = self.boundary_vertices()
        v, e = self.vertex_ring_edges_iterators(sort=True)
        v, vj = self.vertex_ring_vertices_iterators(sort=True)
        cotan = self.edge_cotangents_weigths()
        cotan = np.array([cotan, cotan, cotan]).T
        M = cotan[e] * (self.vertices[vj] - self.vertices[v])
        M = utilities.sum_repeated(M, v)
        M[b,:] = np.array([0,0,0])
        return M

    # -------------------------------------------------------------------------
    #                                  Cutting
    # -------------------------------------------------------------------------

    def cut(self, vertex_index):
        boundary = self.boundary_vertices()
        hc1 = []
        hc2 = []
        if not np.in1d(vertex_index, boundary)[0]:
            return False
        V = self.V
        E = self.E
        W = 2*E
        H = self.halfedges
        v, vj, l = self.vertex_ring_vertices_iterators(return_lengths=True)
        v0 = vertex_index
        vvv = [v0]
        if l[v0] < 3:
            return False
        h = np.where(H[:,0] == v0)[0]
        h = h[np.where(H[h,1] == -1)[0]]
        h0 = copy.copy(h)
        H[h0,0] = V
        n_v1 = np.copy(self.vertices[v0])
        for i in range(l[v0]//2):
            h = H[H[h,4],2]
            H[h,0] = V
        hc1.append(h); hc2.append(H[h,4])
        v1 = H[H[h,4],0][0]
        n_h1 = np.array([[V+1, -1, h0, W+2, h, E]], 'i')
        n_h2 = np.array([[v0, -1, W+3, H[h0,3], H[h,4], H[h,5]]], 'i')
        H[H[h,4],4] = W+1
        H[h,0] = V
        H[h,4] = W
        H[h,5] = E
        H[H[h0,3],2] = W+1
        H[h0,3] = W
        H = np.vstack((H, n_h1, n_h2))
        self._vertices = np.vstack((self.vertices, n_v1))
        W += 2
        E += 1
        V += 1
        while not np.in1d(v1, boundary)[0]:
            v0 = v1
            h = H[h,2]
            H[h,0] = V
            n_v1 = np.copy(self.vertices[v0])
            for i in range(int(l[v0]//2 - 1)):
                h = H[H[h,4],2]
                H[h,0] = V
            hc1.append(h); hc2.append(H[h,4])
            v1 = H[H[h,4],0][0]
            vvv.append(v0)
            n_h1 = np.array([[V+1, -1, W-2, W+2, h, E]], 'i')
            n_h2 = np.array([[v0, -1, W+3, W-1, H[h,4], H[h,5]]], 'i')
            H[H[h,4],4] = W+1
            H[h,0] = V
            H[h,4] = W
            H[h,5] = E
            H = np.vstack((H, n_h1, n_h2))
            self._vertices = np.vstack((self.vertices, n_v1))
            W += 2
            E += 1
            V += 1
        W -= 2
        n_v1 = np.copy(self.vertices[v0])
        while H[h,1] != -1:
            h = H[H[h,2],4]
            H[H[h,4],0] = V
        H[W+1,2] = copy.copy(H[h,2])
        H[H[h,2],3] = W+1
        H[h,2] = W
        H[W,3] = h
        self._halfedges = H
        n_v1 = np.copy(self.vertices[v1])
        self._vertices = np.vstack((self.vertices, n_v1))
        self.topology_update()
        Hcouple = np.array(hc2)
        '''
        from geometrylab import vtkplot
        pl_m = vtkplot.Edges(self)
        pl_p = vtkplot.Points(self.vertices[vvv])
        vtkplot.view([pl_m,pl_p])
        vtkplot.check(self)
        '''
        return Hcouple

    # -------------------------------------------------------------------------
    #                          Edit global connectivity
    # -------------------------------------------------------------------------

    def is_triangular_mesh(self):
        l = self.face_lengths()
        if len(np.where(l != 3)[0] > 0):
            return False
        else:
            return True

    def loop(self, steps=1):
        if not self.is_triangular_mesh():
            return False
        for i in range(steps):
            self._loop()

    def catmull_clark(self, steps=1):
        for i in range(steps):
            self._catmull_clark()

    def _loop(self):
        V = self.V
        H = self.halfedges
        _, h1 = np.unique(H[:,1], True)
        h1 = np.delete(h1, np.where(H[h1,1] == -1))
        h2 = H[h1,2]
        h3 = H[h1,3]
        F0 = np.array((H[h1,5]+V, H[h2,5]+V, H[h3,5]+V)).T
        F1 = np.array((H[h1,0], H[h1,5]+V, H[H[h1,3],5]+V)).T
        F2 = np.array((H[h2,0], H[h2,5]+V, H[H[h2,3],5]+V)).T
        F3 = np.array((H[h3,0], H[h3,5]+V, H[H[h3,3],5]+V)).T
        new_faces = np.vstack((F0, F1, F2, F3)).tolist()
        v, vj, l = self.vertex_ring_vertices_iterators(sort=True,
                                                    return_lengths=True)
        c = 1./l * (5./8 - (3./8 + 1./4*np.cos(2*np.pi*l**(-1.)))**2)
        d = 1 - l*c
        c = np.array([c,c,c]).T
        d = np.array([d,d,d]).T
        ring = utilities.sum_repeated(self.vertices[vj], v)
        V0 = c*ring + d*self.vertices
        _, e = np.unique(H[:,5], True)
        v1 = self.vertices[H[e,0]]
        v2 = self.vertices[H[H[e,4],0]]
        v3 = self.vertices[H[H[e,3],0]]
        v4 = self.vertices[H[H[H[e,4],3],0]]
        be = self.boundary_edges()
        v3[be] = v1[be]
        v4[be] = v2[be]
        V1 = 3./8*v1 + 3./8*v2 + 1./8*v3 + 1./8*v4
        bh = np.where(H[:,1] == -1)[0]
        v0 = self.vertices[H[bh,0]]
        v5 = self.vertices[H[H[bh,3],0]]
        v6 = self.vertices[H[H[bh,2],0]]
        V0[H[bh,0]] = 1./8*v6 + 1./8*v5 + 3./4*v0
        new_vertices = np.vstack((V0, V1))
        self.make_mesh(new_vertices, new_faces)

    def _catmull_clark(self):
        V = self.V
        E = self.E
        H = self.halfedges
        h = self.inner_halfedges()
        v1 = H[h,0]
        v2 = H[h,5] + V
        v3 = H[h,1] + V + E
        v4 = H[H[h,3],5] + V
        new_faces = (np.array([v1, v2, v3, v4]).T).tolist()
        v  = H[:,0]
        i  = np.argsort(v)
        v  = v[i]
        v1 = H[H[:,4],0][i]
        v2 = H[H[H[:,2],4],0][i]
        s = np.ones(len(v))
        l = utilities.sum_repeated(s,v)
        d = 1. / (4*l**2)
        c = 3. / (2*l**2)
        e = 1 - 7./(4*l)
        c = np.array([c,c,c]).T
        d = np.array([d,d,d]).T
        e = np.array([e,e,e]).T
        v1 = utilities.sum_repeated(self.vertices[v1], v)
        v2 = utilities.sum_repeated(self.vertices[v2], v)
        old_vertices = c*v1 + d*v2 + e*self.vertices
        _, e = np.unique(H[:,5], True)
        v1 = self.vertices[H[e,0]]
        v2 = self.vertices[H[H[e,4],0]]
        v3 = self.vertices[H[H[e,3],0]]
        v4 = self.vertices[H[H[H[e,4],3],0]]
        v5 = self.vertices[H[H[H[e,2],2],0]]
        v6 = self.vertices[H[H[H[H[e,4],2],2],0]]
        be = self.boundary_edges()
        v3[be] = v5[be] = v1[be]
        v4[be] = v6[be] = v2[be]
        mid_points = 3./8*(v1 + v2) + 1./16*(v3 + v4 + v5 + v6)
        bh = np.where(H[:,1] == -1)[0]
        v0 = self.vertices[H[bh,0]]
        v5 = self.vertices[H[H[bh,3],0]]
        v6 = self.vertices[H[H[bh,2],0]]
        old_vertices[H[bh,0]] = 1./8*v6 + 1./8*v5 + 3./4*v0
        barycenters = self.face_barycenters()
        new_vertices = np.vstack((old_vertices, mid_points, barycenters))
        self.make_mesh(new_vertices, new_faces)

    def dual_mesh(self, make_boundary=True):
        H = self.halfedges
        HD = np.copy(H)
        b = np.where(H[:,1] == -1)[0]
        B = b.shape[0]
        fb = np.arange(B) + self.F
        hb1 = np.arange(B) + 2*self.E
        hb2 = np.arange(B) + B + 2*self.E
        eb = np.arange(B) + self.E
        HD[b,1] = np.copy(fb)
        HD[b,3] = np.copy(hb1)
        HD[H[b,3],2] = np.copy(hb2)
        Hb1 = np.zeros((B,6), dtype=int)
        Hb2 = np.zeros((B,6), dtype=int)
        Hb1[:,0] = -1
        Hb2[:,0] = H[b,0]
        Hb1[:,1] = fb
        Hb2[:,1] = HD[H[b,3],1]
        Hb1[:,2] = b
        Hb2[:,2] = HD[H[b,3],3]#HD[H[b,2],2]##
        Hb1[:,3] = HD[b,2]
        Hb2[:,3] = H[b,3]
        Hb1[:,4] = hb2
        Hb2[:,4] = hb1
        Hb1[:,5] = eb
        Hb2[:,5] = eb
        HD = np.vstack((HD, Hb1, Hb2))
        HR = np.copy(HD)
        HD[:,0] = HR[HR[:,4],1]
        HD[:,1] = HR[:,0]
        HD[:,2] = HR[HR[:,3],4]
        HD[:,3] = HR[HR[:,4],2]
        dual_vertices = self.face_barycenters()
        face_normals = self.face_normals()
        edge_vec = self.vertices[H[b,0]] - self.vertices[H[H[b,2],0]]
        normals = np.cross(edge_vec, face_normals[H[H[b,4],1]])
        new_vertices = self.edge_mid_points()[H[b,5]] + normals/2
        dual_vertices = np.vstack((dual_vertices, new_vertices))
        dual_mesh = Mesh()
        dual_mesh._halfedges = HD
        dual_mesh._vertices = dual_vertices
        dual_mesh.topology_update()
        return dual_mesh

    def delete_faces(self, faces):
        if type(faces) is int:
            faces = [faces]
        H = self.halfedges
        self._open_trash()
        self._cull_faces(faces)
        hf = np.arange(H.shape[0])[np.in1d(H[:,1], faces)]
        bh = hf[H[H[hf,4],1] == -1]
        self._cull_halfedges(bh)
        self._cull_halfedges(H[bh,4])
        self._cull_edges(H[bh,5])
        dh = hf[np.in1d(H[hf,4], hf)]
        self._cull_halfedges(dh)
        self._cull_halfedges(H[dh,4])
        self._cull_edges(H[dh,5])
        H[hf,1] = -1
        self._clean()
        self.delete_unconnected_vertices()

    def delete_unconnected_vertices(self):
        H = self.halfedges
        v = np.arange(len(self.vertices))
        cull = np.invert(np.in1d(v, H[:,0]))
        self._open_trash()
        self._cull_vertices(v[cull])
        self._clean()

    def exploded_mesh(self):
        f, v = self.face_vertices_iterators()
        vertices = self.vertices[v]
        k = np.arange(len(v))
        faces_list = [[] for i in range(self.F)]
        for i in range(len(f)):
            faces_list[f[i]].append(k[i])
        exp_mesh = Mesh()
        exp_mesh.make_mesh(vertices, faces_list)
        return exp_mesh

    # -------------------------------------------------------------------------
    #                            Edit edge connectivity
    # -------------------------------------------------------------------------

    def delete_edge(self, edge_index):
        h = self.edge_halfedge(edge_index)
        self._open_trash()
        self._delete_halfedge(h)
        self._clean()

    def flip_edge(self, edge_index):
        h = self.edge_halfedge(edge_index)
        if not self.is_halfedge_bounding_tri_faces(h):
            return False
        self._flip_halfedge(h)
        self.topology_update()

    def split_edge(self, edge_index):
        h = self.edge_halfedge(edge_index)
        if not self.is_halfedge_bounding_tri_faces(h):
            return False
        self._expand_arrays()
        self._split_halfedge(h)
        self._cull_arrays()

    def collapse_edge(self, edge_index):
        h = self.edge_halfedge(edge_index)
        if not self.is_halfedge_bounding_tri_faces(h):
            return False
        self._open_trash()
        self._collapse_halfedge(h)
        self._clean()

    # -------------------------------------------------------------------------
    #                                   Remesh
    # -------------------------------------------------------------------------

    def equalize_valences(self):
        if not self.is_triangular_mesh():
            return False
        H = self.halfedges
        _, he = np.unique(H[:,5], return_index=True)
        t = np.repeat(6, self.V)
        t[self.boundary_vertices()] = 4
        t[self.mesh_corners()] = 3
        _, _, l = self.vertex_ring_vertices_iterators(True, False, True)
        for h in he:
            a = H[h,0]
            b = H[H[h,4],0]
            c = H[H[h,3],0]
            d = H[H[H[h,4],3],0]
            deviation_0  = (l[a] - t[H[h,0]])**2
            deviation_0 += (l[b] - t[H[H[h,4],0]])**2
            deviation_0 += (l[c] - t[H[H[h,3],0]])**2
            deviation_0 += (l[d] - t[H[H[H[h,4],3],0]])**2
            deviation_1  = (l[a] - t[H[h,0]] - 1)**2
            deviation_1 += (l[b] - t[H[H[h,4],0]] - 1)**2
            deviation_1 += (l[c] - t[H[H[h,3],0]] + 1)**2
            deviation_1 += (l[d] - t[H[H[H[h,4],3],0]] + 1)**2
            if deviation_1 < deviation_0:
                if self._flip_halfedge(h):
                    l[a] -= 1
                    l[b] -= 1
                    l[c] += 1
                    l[d] += 1
        self.topology_update()
        return True

    def split_edges(self, max_length):
        if not self.is_triangular_mesh():
            return False
        H = self.halfedges
        _, he = np.unique(H[:,5], return_index=True)
        self._expand_arrays(len(he))
        for h in he:
            if self.halfedge_length(h) > max_length:
                self._split_halfedge(h)
        self._cull_arrays()
        return True

    def collapse_edges(self, min_length):
        if not self.is_triangular_mesh():
            return False
        H = self.halfedges
        _, he = np.unique(H[:,5], return_index=True)
        self._open_trash()
        for h in he:
            if self._is_culled_halfedge(h):
                h = H[h,4]
            if self.halfedge_length(h) < min_length:
                self._collapse_halfedge(h)
        self._clean()
        return True

    # -------------------------------------------------------------------------
    #                             Local connectivity
    # -------------------------------------------------------------------------

    def edge_halfedge(self, edge_index):
        H = self.halfedges
        h = np.where(H[:,5] == edge_index)[0][0]
        return h

    def vertex_halfedge(self, vertex_index):
        H = self.halfedges
        v = np.where(H[:,0] == vertex_index)[0][0]
        return v

    def halfedge_ring(self, halfedge_index):
        H = self.halfedges
        h0 = halfedge_index
        ring = [h0]
        h = H[H[h0,3],4]
        while h != h0:
            ring.append(h)
            h = H[H[h,3],4]
        return ring

    def vertex_ring_vertices(self, vertex_index):
        h = self.vertex_halfedge(vertex_index)
        ring = self.halfedge_ring_vertices(h)
        return ring

    def vertex_multiple_ring_vertices(self, vertex_index, depth=1):
        vi, vj = self.vertex_ring_vertices_iterators()
        ring = np.array([], dtype='i')
        search = np.array([vertex_index], dtype='i')
        for i in range(int(depth)):
            vring = np.array([], dtype='i')
            for v in search:
                vring = np.hstack((vj[vi == v], vring))
            vring = np.unique(vring)
            vring = vring[np.invert(np.in1d(vring, ring))]
            search = vring
            ring = np.hstack((ring, vring))
            if len(ring) == self.V:
                return ring
        return np.unique(ring)

    def halfedge_ring_vertices(self, halfedge_index):
        H = self.halfedges
        ring = self.halfedge_ring(halfedge_index)
        vertices = H[H[ring,2],0]
        return vertices

    def halfedge_ring_faces(self, halfedge_index):
        H = self.halfedges
        ring = self.halfedge_ring(halfedge_index)
        faces = H[H[ring,2],1]
        return faces

    def halfedge_face_vertices(self, halfedge_index):
        H = self.halfedges
        ring = self.halfedge_face_ring(halfedge_index)
        vertices = H[ring,0]
        return vertices

    # -------------------------------------------------------------------------
    #                             Local queries
    # -------------------------------------------------------------------------

    def is_boundary_halfedge_ring(self, ring):
        H = self.halfedges
        for h in ring:
            if H[h,1] == -1:
                v0 = H[h,0]
                v1 = H[H[h,2],0]
                v2 = H[H[h,3],0]
                E1 = (self.vertices[v1] - self.vertices[v0])
                E2 = (self.vertices[v0] - self.vertices[v2])
                E1 = E1 / (E1[0]**2 + E1[1]**2 + E1[2]**2 + 1e-10)**0.5
                E2 = E2 / (E2[0]**2 + E2[1]**2 + E2[2]**2 + 1e-10)**0.5
                dot = E1[0]*E2[0] + E1[1]*E2[1] + E1[2]*E2[2]
                if dot < self.corner_tolerance:
                    return 2
                else:
                    return 1
        return 0

    def is_halfedge_bounding_tri_faces(self, halfedge_index):
        H = self.halfedges
        h = halfedge_index
        for i in range(2):
            counter = 1
            h0 = h
            h = H[h,2]
            while h != h0:
                h = H[h,2]
                counter += 1
                if counter > 3:
                    return False
            h = H[halfedge_index,4]
        return True

    # -------------------------------------------------------------------------
    #                             Local geometry
    # -------------------------------------------------------------------------

    def halfedge_length(self, halfedge_index):
        H = self.halfedges
        h = halfedge_index
        V1 = self.vertices[H[h,0]]
        V2 = self.vertices[H[H[h,4],0]]
        E = V1 - V2
        return (E[0]**2 + E[1]**2 + E[2]**2)**0.5

    # -------------------------------------------------------------------------
    #                         Edit half-edge connectivity
    # -------------------------------------------------------------------------

    def _delete_halfedge(self, halfedge_index):
        H = self.halfedges
        h1 = int(halfedge_index)
        h2 = H[h1,4]
        f1 = H[h1,1]
        f2 = H[h2,1]
        if f1 == -1 or f2 == -1:
            return False
        h3 = H[h1,2]
        h4 = H[h1,3]
        h5 = H[h2,2]
        h6 = H[h2,3]
        H[h3,3] = h6
        H[h6,2] = h3
        H[h5,3] = h4
        H[h4,2] = h5
        H[np.where(H[:,1] == f2)[0],1] = f1
        self._cull_halfedge(h1)
        self._cull_halfedge(h2)
        self._cull_edge(H[h1,5])
        self._cull_face(f2)
        self._F -= 1
        self._E -= 1
        return True

    def _flip_halfedge(self, halfedge_index):
        H = self.halfedges
        h1 = int(halfedge_index)
        h2 = H[h1,4]
        f1 = H[h1,1]
        f2 = H[h2,1]
        if f1 == -1 or f2 == -1:
            return False
        ring1 = self.halfedge_ring_vertices(h1)
        ring2 = self.halfedge_ring_vertices(h2)
        if len(ring1) <= 3 or len(ring2) <= 3:
            return False
        v3 = H[H[h1,3],0]
        v4 = H[H[h2,3],0]
        ring3 = self.halfedge_ring_vertices(H[h1,3])
        for h in ring3:
            vertices = [H[h,0], H[H[h,4],0]]
            if v3 in vertices and v4 in vertices: #?
                return False
        h3 = H[h1,2]
        h4 = H[h1,3]
        h5 = H[h2,2]
        h6 = H[h2,3]
        i = [h1,h1,h1,h2,h2,h2,h3,h3,h4,h4,h5,h5,h6,h6,h6,h4]
        j = [0, 3, 2, 0, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1, 1 ]
        d = [v3,h3,h6,v4,h5,h4,h6,h1,h2,h5,h4,h2,h1,h3,f1,f2]
        H[i,j] = np.array(d)
        return True

    def _split_halfedge(self, halfedge_index):
        V = self.V
        E = self.E
        F = self.F
        K = int(2*E)
        H = self.halfedges
        h1 = int(halfedge_index)
        h2 = H[h1,4]
        v1 = H[h1,0]
        v2 = H[h2,0]
        f1 = H[h1,1]
        f2 = H[h2,1]
        H1 = np.copy(H[h1])
        H2 = np.copy(H[h2])
        new_vertex = (self.vertices[v1] + self.vertices[v2]) / 2
        H[h2,0] = V
        H[h1,2] = K
        H[h2,3] = K+1
        H[H2[3],2] = K+1
        H[H1[2],3] = K
        nh1 = K + 0
        nh2 = K + 1
        H[K,:] = [V, f1, H1[2], h1, K+1, E]
        H[K+1,:] = [v2, f2 , h2, H2[3], K, E]
        K0 = K + 0
        K += 2
        E += 1
        if f1 != -1:
            v3 = H[H1[3],0]
            H[K0,1] = F
            H[K0,3] = K+1
            H[H[K0,2],1] = F
            H[H[K0,2],2] = K+1
            H[H[K0,2],3] = nh1
            H[H1[3],3] = K
            H[h1,2] = K
            H[K,:] = [V, f1, H1[3], h1, K+1, E]
            H[K+1,:] = [v3, F, nh1, H1[2], K, E]
            K += 2
            F += 1
            E += 1
        if f2 != -1:
            v4 = H[H2[3],0]
            H[K0+1,1] = F
            H[K0+1,2] = K+1
            H[H[K0+1,3],1] = F
            H[H[K0+1,3],3] = K+1
            H[H[K0+1,3],2] = nh2
            H[H2[2],2] = K
            H[h2,3] = K
            H[K,:] = [v4, f2, h2, H2[2], K+1, E]
            H[K+1,:] = [V, F, H2[3], nh2, K, E]
            F += 1
            E += 1
        self._halfedges = H
        self._vertices[V,:] = new_vertex
        self._F = F
        self._E = E
        self._V += 1
        return True

    def _collapse_halfedge(self, halfedge_index):
        H = self.halfedges
        h1 = halfedge_index
        if self._is_culled_halfedge(h1):
            return False
        e = H[h1, 5]
        h2 = H[h1,4]
        H1 = np.copy(H[h1])
        H2 = np.copy(H[h2])
        h3 = H1[2]
        h4 = H1[3]
        h5 = H2[2]
        h6 = H2[3]
        h7 = H[h3,4]
        h8 = H[h4,4]
        h9 = H[h5,4]
        h10 = H[h6,4]
        f1 = H1[1]
        f2 = H2[1]
        v1 = H[h1,0]
        v2 = H[h2,0]
        ring1 = self.halfedge_ring(h1)
        ring2 = self.halfedge_ring(h2)
        vring1 = H[H[ring1,2],0]
        vring2 = H[H[ring2,2],0]
        counter = 0
        for v in vring1:
            if v in vring2:
                counter += 1
        if f1 != -1 and f2 != -1 and counter > 2:
            return False
        if f1 == -1 and f2 != -1 and counter > 1:
            return False
        if f2 == -1 and f1 != -1 and counter > 1:
            return False
        V1 = np.copy(self.vertices[v1])
        V2 = np.copy(self.vertices[v2])
        bv1 = self.is_boundary_halfedge_ring(ring1)
        bv2 = self.is_boundary_halfedge_ring(ring2)
        if bv1:
            if not bv2:
                new_v = V1
            elif f1 != -1 and f2 != -1:
                return False
            elif bv1 == 2:
                new_v = V1
            elif bv2 == 2:
                new_v = V2
            else:
                new_v = (V1 + V2) / 2
        elif bv2:
            new_v = V2
        else:
           new_v = (V1 + V2) / 2
        H[ring2,0] = v1
        self._vertices[v1] = new_v
        self._cull_halfedge(h1)
        self._cull_halfedge(h2)
        self._cull_edge(e)
        self._cull_vertex(v2)
        self._V -= 1
        self._E -= 1
        if f1 != -1:
            H[h7,4] = h8
            H[h8,4] = h7
            H[h7,5] = H[h8,5]
            self._cull_halfedge(h3)
            self._cull_halfedge(h4)
            self._cull_edge(H[h3,5])
            self._cull_face(f1)
            self._F -= 1
            self._E -= 1
        else:
            H[h4,2] = h3
            H[h3,3] = h4
        if f2 != -1:
            H[h9,4] = h10
            H[h10,4] = h9
            H[h10,5] = H[h9,5]
            self._cull_halfedge(h5)
            self._cull_halfedge(h6)
            self._cull_edge(H[h6,5])
            self._cull_face(f2)
            self._F -= 1
            self._E -= 1
        else:
            H[h6,2] = h5
            H[h5,3] = h6
        return True

    #--------------------------------------------------------------------------
    #                          Culling and Reallocation
    #--------------------------------------------------------------------------

    def _is_culled_halfedge(self, halfedge_index):
        return self.__Htrash[halfedge_index]

    def _cull_halfedge(self, halfedge_index):
        self.__Htrash[halfedge_index] = True

    def _cull_face(self, face_index):
        self.__Ftrash.append(face_index)

    def _cull_edge(self, edge_index):
        self.__Etrash.append(edge_index)

    def _cull_vertex(self, vertex_index):
        self.__Vtrash.append(vertex_index)

    def _cull_halfedges(self, halfedge_indices):
        self._cull_halfedge(halfedge_indices)

    def _cull_faces(self, face_indices):
        try:
            face_indices = face_indices.tolist()
        except:
            pass
        self.__Ftrash.extend(face_indices)

    def _cull_edges(self, edge_indices):
        try:
            edge_indices = edge_indices.tolist()
        except:
            pass
        self.__Etrash.extend(edge_indices)

    def _cull_vertices(self, vertex_indices):
        try:
            pass#vertex_indices = vertex_indices.tolist()
        except:
            pass
        self.__Vtrash.extend(vertex_indices)

    def _expand_arrays(self, n=1):
        NH = -np.ones((6*n, 6), 'i')
        NV = np.empty((n, 3))
        NV.fill(np.nan)
        self._halfedges = np.vstack((self._halfedges, NH))
        self._vertices = np.vstack((self._vertices, NV))

    def _cull_arrays(self):
        H = self.halfedges
        Hdelete = np.where(H[:,0] == -1)[0]
        self._halfedges = np.delete(H, Hdelete, axis=0)
        V = self.vertices
        self._vertices = np.delete(V, np.where(np.isnan(V[:,0]))[0], axis=0)
        self.topology_update()

    def _open_trash(self):
        self.__Htrash = np.zeros(self.halfedges.shape[0], 'b')
        self.__Vtrash = []
        self.__Ftrash = []
        self.__Etrash = []

    def _close_trash(self):
        self.__Htrash = None
        self.__Vtrash = None
        self.__Ftrash = None
        self.__Etrash = None

    def _clean(self):
        H = self.halfedges
        if self.__Htrash is None:
            return False
        self._V = np.amax(self.halfedges[:,0] + 1)
        self._F = np.amax(self.halfedges[:,1] + 1)
        self._E = np.amax(self.halfedges[:,5] + 1)
        self.__Htrash = np.nonzero(self.__Htrash)[0]
        self.__Htrash = np.unique(self.__Htrash)
        self.__Vtrash = np.unique(self.__Vtrash)
        self.__Ftrash = np.unique(self.__Ftrash)
        self.__Etrash = np.unique(self.__Etrash)
        self.__Vtrash = sorted(self.__Vtrash)
        self.__Ftrash = sorted(self.__Ftrash)
        self.__Etrash = sorted(self.__Etrash)
        try:
            V = max(self.V, max(self.__Vtrash))
        except:
            V = self.V
        try:
            F = max(self.F, max(self.__Ftrash))
        except:
            F = self.F
        try:
            E = max(self.E, max(self.__Etrash))
        except:
            E = self.E
        if len(self.__Htrash) > 0:
            roll = np.roll(self.__Htrash, 1)
            roll[0] = 0
            ind = np.arange(len(self.__Htrash))
            hd = np.repeat(ind, self.__Htrash - roll)
            hd = np.hstack((hd, np.repeat(hd[-1]+1, 2*E - self.__Htrash[-1])))
            h = np.arange(2*E) - hd
            H[:,2] = h[H[:,2]]
            H[:,3] = h[H[:,3]]
            H[:,4] = h[H[:,4]]
            H = np.delete(H, self.__Htrash, axis=0)
        if len(self.__Vtrash) > 0:
            roll = np.roll(self.__Vtrash, 1)
            roll[0] = 0
            ind = np.arange(len(self.__Vtrash))
            vd = np.repeat(ind, self.__Vtrash - roll)
            vd = np.hstack((vd, np.repeat(vd[-1] + 1, V - self.__Vtrash[-1])))
            v = np.arange(V) - vd
            H[:,0] = v[H[:,0]]
        if len(self.__Ftrash) > 0:
            roll = np.roll(self.__Ftrash, 1)
            roll[0] = 0
            ind = np.arange(len(self.__Ftrash))
            fd = np.repeat(ind, self.__Ftrash - roll)
            fd = np.hstack((fd, np.repeat(fd[-1] + 1, F - self.__Ftrash[-1])))
            f = np.arange(F) - fd
            f = np.hstack((f, [-1]))
            H[:,1] = f[H[:,1]]
        if len(self.__Etrash) > 0:
            roll = np.roll(self.__Etrash, 1)
            roll[0] = 0
            ind = np.arange(len(self.__Etrash))
            ed = np.repeat(ind, self.__Etrash - roll)
            ed = np.hstack((ed, np.repeat(ed[-1] + 1, E - self.__Etrash[-1])))
            e = np.arange(E) - ed
            H[:,5] = e[H[:,5]]
        self._vertices = np.delete(self.vertices, self.__Vtrash, axis=0)
        self._halfedges = H
        self._close_trash()
        self.topology_update()
        return True

    def _connectivity_check(self):
        H = self.halfedges
        out = '*** Mesh connectivity check ***\n'
        for h in range(len(H)):
            if H[H[h,4],4] != h:
                out += ('halfedge {}: twin error\n').format(h)
            if H[H[h,4],5] != H[h,5]:
                out += ('halfedge {}: edge error\n').format(h)
            if H[H[h,2],1] != H[h,1]:
                out += ('halfedge {}: next face error\n').format(h)
            if H[H[h,3],1] != H[h,1]:
                out += ('halfedge {}: previous face error\n').format(h)
            if H[H[h,2],0] != H[H[h,4],0]:
                out += ('halfedge {}: next twin origin error\n').format(h)
            if H[H[H[h,3],4],0] != H[h,0]:
                out += ('halfedge {}: previous twin origin error\n').format(h)
        for v in range(self.V):
            origins = np.where(H[:,0] == v)[0]
            if len(origins) < 2:
                out += ('vertex {}: incoming edges < 2\n').format(v)
        for e in range(self.E):
            halfedges = np.where(H[:,5] == e)[0]
            if len(halfedges) != 2:
                out += ('edge {}: connectivity error\n').format(v)
        out += 'check finished'



