#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from traits.api import HasTraits, Instance, Property, Enum, Button,String,\
                       on_trait_change, Float, Bool, Int,Constant, ReadOnly,\
                       List, Array

from traitsui.api import View, Item, HSplit, VSplit, InstanceEditor, HGroup,\
                         Group, ListEditor, Tabbed, VGroup, CheckListEditor,\
                         ArrayEditor, Action, ToolBar

from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel

from mayavi.core.ui.mayavi_scene import MayaviScene

from pyface.image_resource import ImageResource

#------------------------------------------------------------------------------

from geometrylab.gui.scenemanager import SceneManager

from geometrylab.vtkplot.pointsource import Points

from geometrylab.vtkplot.vectorsource import Vectors

from geometrylab.vtkplot.polylinesource import Polyline

# -----------------------------------------------------------------------------

'''check.py: an interactive checker for mesh connectivity and normals'''

__author__ = 'Davide Pellis'


###############################################################################
#################################### VIEWER ###################################
###############################################################################

class MeshCheck(HasTraits):

    scene = Instance(MlabSceneModel, ())

    editor = SceneEditor(scene_class=MayaviScene)

    vertex_normals = Bool(False)

    face_normals = Bool(False)

    edge_normals = Bool(False)

    boundary_normals = Bool(False)

    half_edge = Bool(False)

    vertex_ring = Bool(False)

    faces = Bool(False)

    flip_edges = Bool(False)

    collapse_edges = Bool(False)

    split_edges = Bool(False)

    boundary_polyline = Bool(False)

    move_vertices = Bool(False)

    select_faces = Bool(False)

    select = Bool(False)

    multiple_ring = Bool(False)

    depth = Int(1)

    #--------------------------------------------------------------------------
    view = View(
                HSplit(VGroup(
                              Item('faces'),
                              Item('vertex_normals'),
                              Item('face_normals'),
                              Item('edge_normals'),
                              Item('boundary_normals'),
                              Item('half_edge'),
                              Item('vertex_ring'),
                              Item('flip_edges'),
                              Item('collapse_edges'),
                              Item('split_edges'),
                              Item('boundary_polyline'),
                              Item('move_vertices'),
                              Item('select_faces'),
                              Item('select'),
                              Item('multiple_ring'),
                              Item('depth'),
                              show_border=True),
                       Item('scene',
                            editor=editor,
                            show_label=False,
                            resizable=True,
                            height=800,
                            width=1200,
                            ),
                        ),
                resizable=True,
                title = 'Mesh Check',
                icon = ImageResource('logo.png'),
                )
    #--------------------------------------------------------------------------

    def __init__(self, mesh):
        HasTraits.__init__(self)
        self.scenemanager = SceneManager(scene_model=self.scene)
        self.scenemanager.add_object(mesh)
        self.mesh_manager.plot_edges()

    @property
    def mesh(self):
        return self.mesh_manager.mesh

    @property
    def mesh_manager(self):
        return self.scenemanager.current_object

    def start(self):
        self.configure_traits()

    def set_state(self, name):
        if name != 'half_edge':
            self.half_edge = False
        if name != 'vertex_ring':
            self.vertex_ring = False
        if name != 'multiple_ring':
            self.multiple_ring = False
        if name != 'flip_edges':
            self.flip_edges = False
        if name != 'collapse_edges':
            self.collapse_edges = False
        if name != 'split_edges':
            self.split_edges = False
        if name != 'move_vertices':
            self.move_vertices = False
        if name != 'select_faces':
            self.select_faces = False
        if name != 'select':
            self.select = False

    #--------------------------------------------------------------------------

    @on_trait_change('select')
    def select_mesh(self):
        if self.select:
            self.set_state('select')
            self.scenemanager.select_object()
        else:
            self.scenemanager.select_off()

    @on_trait_change('move_vertices')
    def move_vertices_interactive(self):
        if self.move_vertices:
            self.set_state('move_vertices')
            def callback():
                self.mesh_manager.plot_edges()
                self.mesh_manager.plot_faces()
            self.mesh_manager.move_vertices(callback)
        else:
            self.mesh_manager.move_vertices_off()

    @on_trait_change('vertex_normals')
    def plot_vertex_normals(self, mesh):
        if self.vertex_normals:
            scale = 0.3
            Vn = Vectors(self.mesh.vertex_normals(),
                           self.mesh,
                           anchor_mode='vertex',
                           scale_factor = scale,
                           position = 'tail',
                           color = 'b',
                           name = 'vn')
            self.mesh_manager.add([Vn])
        else:
            self.mesh_manager.remove('vn')

    @on_trait_change('face_normals')
    def plot_face_normals(self, mesh):
        if self.face_normals:
            scale = 0.3
            Vf = Vectors(self.mesh.face_normals(),
                           self.mesh,
                           anchor_mode='face',
                           scale_factor = scale,
                           position = 'tail',
                           color = 'r',
                           name = 'vf')
            self.mesh_manager.add([Vf])
        else:
            self.mesh_manager.remove('vf')

    @on_trait_change('edge_normals')
    def plot_edge_normals(self, mesh):
        if self.edge_normals:
            scale = 0.3
            Ve = Vectors(self.mesh.edge_normals(),
                           self.mesh,
                           anchor_mode='edge',
                           scale_factor = scale,
                           position = 'tail',
                           color = 'g',
                           name = 've')
            self.mesh_manager.add([Ve])
        else:
            self.mesh_manager.remove('ve')

    @on_trait_change('boundary_normals')
    def plot_boundary_normals(self, mesh):
        if self.boundary_normals:
            scale = 0.3
            Vb = Vectors(self.mesh.boundary_normals(),
                           self.mesh,
                           anchor_mode='vertex',
                           scale_factor = scale,
                           position = 'tail',
                           color = 'm',
                           name = 'vb')
            self.mesh_manager.add([Vb])
        else:
            self.mesh_manager.remove('vb')

    def plot_ring(self, vertex_index):
        v, ring = self.mesh.vertex_ring_vertices_iterators(order=True)
        ring = ring[np.where(v == vertex_index)[0]]
        data = np.arange(len(ring))
        P = Points(self.mesh.vertices[ring],
                   point_data = data,
                   color = 'bwr',
                   lut_range = '-:+',
                   radius = self.mesh_manager.r*3,
                   name = 'points')
        self.mesh_manager.add(P)

    def plot_multiple_ring(self, vertex_index):
        ring = self.mesh.vertex_ring_expansion(vertex_index, depth=self.depth)
        data = np.arange(len(ring))
        P = Points(self.mesh.vertices[ring],
                   point_data = data,
                   color = 'r',
                   radius = self.mesh_manager.r*3,
                   name = 'points')
        self.mesh_manager.add(P)

    def plot_half_edge(self, edge_index):
        H = self.mesh.halfedges
        e = np.where(H[:,5] == edge_index)[0]
        h1 = e[0]; h2 = e[1]
        t1 = H[H[h1,4],5]
        t2 = H[H[h2,4],5]
        edges = np.unique(np.array([edge_index, t1, t2]))
        f1 = H[h1,1]
        f2 = H[h2,1]
        if f1 == -1:
            if f2 != -1:
                faces = np.array([f2])
        elif f2 == -1:
            if f1 != -1:
                faces = np.array([f1])
        else:
            faces = np.array([f1,f2])
        nex1 = H[H[h1,2],5]
        nex2 = H[H[h2,2],5]
        nex = np.array([nex1,nex2])
        pre1 = H[H[h1,3],5]
        pre2 = H[H[h2,3],5]
        pre = np.array([pre1,pre2])
        v1 = H[h1,0]
        v2 = H[h2,0]
        vertices = np.array([v1,v2])
        edge_data = np.zeros(self.mesh.E)
        edge_data[edges] = 1
        edge_data[pre] = 2
        edge_data[nex] = 3
        face_data = np.ones(self.mesh.F)
        face_data[faces] = 0
        self.mesh_manager.plot_edges(edge_data = edge_data,
                                   color = [(0.6,0.6,0.6),'yellow','b','r'],
                                   lut_range = '-:+',
                                   )
        self.mesh_manager.plot_faces(face_data = face_data,
                                   color = ['white',(0.6,0.6,0.6)],
                                   lut_range = '-:+')
        P = Points(self.mesh,
                   vertex_indices = vertices,
                   color = 'yellow',
                   name = 'points',
                   radius = self.mesh_manager.r*3)
        self.mesh_manager.add(P)

    @on_trait_change('half_edge')
    def half_edge_check(self):
        if self.half_edge:
            self.set_state('half_edge')
            self.mesh_manager.on_edge_selection(self.plot_half_edge)
        else:
            self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
            self.mesh_manager.plot_faces()
            self.mesh_manager.remove('points')
            self.mesh_manager.select_edges_off()

    @on_trait_change('vertex_ring')
    def vertex_ring_check(self):
        if self.vertex_ring:
            self.set_state('vertex_ring')
            self.mesh_manager.on_vertex_selection(self.plot_ring)
        else:
            self.mesh_manager.remove('points')
            self.mesh_manager.select_vertices_off()

    @on_trait_change('multiple_ring')
    def vertex_multiple_ring_check(self):
        if self.multiple_ring:
            self.set_state('multiple_ring')
            self.mesh_manager.on_vertex_selection(self.plot_multiple_ring)
        else:
            self.mesh_manager.remove('points')
            self.mesh_manager.select_vertices_off()

    @on_trait_change('faces')
    def face_order(self):
        if self.faces:
            data = np.arange(self.mesh.F)
            self.mesh_manager.plot_faces(face_data=data,
                                       color = 'Yl_or_br')
        else:
            self.mesh_manager.plot_faces()

    def flip(self, edge_index):
        self.mesh.flip_edge(edge_index)
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))

    def collapse(self, edge_index):
        self.mesh.collapse_edge(edge_index)
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))

    def split(self, edge_index):
        self.mesh.split_edge(edge_index)
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))

    @on_trait_change('flip_edges')
    def flip_edges(self):
        if self.flip_edges:
            self.set_state('flip_edges')
            self.mesh_manager.on_edge_selection(self.flip)
        else:
            self.mesh_manager.select_edges_off()

    @on_trait_change('collapse_edges')
    def collapse_edges(self):
        if self.collapse_edges:
            self.set_state('collapse_edges')
            self.mesh_manager.on_edge_selection(self.collapse)
        else:
            self.mesh_manager.select_edges_off()

    @on_trait_change('split_edges')
    def split_edges(self):
        if self.split_edges:
            self.set_state('split_edges')
            self.mesh_manager.on_edge_selection(self.split)
        else:
            self.mesh_manager.select_edges_off()

    @on_trait_change('select_faces')
    def select_faces(self):
        if self.select_faces:
            self.set_state('select_faces')
            self.mesh_manager.select_faces()
        else:
            self.mesh_manager.select_faces_off()




# -----------------------------------------------------------------------------

def check(mesh):
    m = MeshCheck(mesh)
    m.start()
