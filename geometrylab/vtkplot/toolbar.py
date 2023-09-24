#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from traits.api import HasTraits, Instance, Property, Enum, Button,String,\
                       on_trait_change, Float, Bool, Int,Constant, ReadOnly,\
                       List, Array, Range

from traitsui.api import View, Item, HSplit, VSplit, InstanceEditor, HGroup,\
                         Group, ListEditor, Tabbed, VGroup, CheckListEditor,\
                         ArrayEditor, Action, ToolBar

from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel

from mayavi.core.ui.mayavi_scene import MayaviScene

from pyface.image_resource import ImageResource

#------------------------------------------------------------------------------

from geometrylab.vtkplot.scenemanager import SceneManager

from geometrylab.vtkplot.pointsource import Points

from geometrylab.vtkplot.vectorsource import Vectors

from geometrylab.vtkplot.polylinesource import Polyline

# -----------------------------------------------------------------------------

'''check.py: an interactive checker for mesh connectivity and normals'''

__author__ = 'Davide Pellis'


###############################################################################
#################################### VIEWER ###################################
###############################################################################

class IncrementalRemesh(HasTraits):

    remesh_factor = Range(0.1, 1.5, 1)

    corner_tolerance = Range(-1.0, 1.0, 0.3)

    remesh_button = Button(label='Remesh')

    view = View(VGroup('remesh_factor',
                       'corner_tolerance',
                       Item('remesh_button',
                            show_label=False),
                       show_border=True),
                title = 'Incremental Remesh',
                icon = ImageResource('logo.png'))

    def __init__(self):
        HasTraits.__init__(self)
        self._mesh_manager = None

    @property
    def mesh(self):
        return self.mesh_manager.mesh

    @property
    def mesh_manager(self):
        return self._mesh_manager

    @mesh_manager.setter
    def mesh_manager(self, mesh_manager):
        self._mesh_manager = mesh_manager
        self.corner_tolerance = self.mesh.corner_tolerance

    def start(self):
        self.configure_traits()

    def update_plot(self):
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))

    @on_trait_change('corner_tolerance')
    def corner_tolerance(self):
        self.mesh.corner_tolerance = self.corner_tolerance

    @on_trait_change('remesh_button')
    def incremental_remesh(self):
        try:
            self.mesh.incremental_remesh(self.remesh_factor)
            self.mesh_manager.update_plot()
        except:
            pass


class MeshEditToolbar(HasTraits):

    scene = Instance(MlabSceneModel, ())

    editor = SceneEditor(scene_class=MayaviScene)

    remesher = IncrementalRemesh()

    flip_edges = Action(action='flip_edges',
                        image=ImageResource('flip.png'),
                        style='toggle',
                        tooltip = 'Flip mesh edges',
                        show_label=False)

    split_edges = Action(action='split_edges',
                         image=ImageResource('split.png'),
                         style='toggle',
                         tooltip = 'Split mesh edges',
                         show_label=False)

    collapse_edges = Action(action='collapse_edges',
                            image=ImageResource('collapse.png'),
                            style='toggle',
                            tooltip = 'Collapse mesh edges',
                            show_label=False)

    catmull_clark = Action(action='catmull_clark',
                            image=ImageResource('catmullclark.png'),
                            style='push',
                            tooltip = 'Catmull-Clark subdivision',
                            show_label=False)

    loop_subdivision = Action(action='loop',
                             image=ImageResource('loop.png'),
                             style='push',
                             tooltip = 'Loop subdivision',
                             show_label=False)

    incremental_remesh = Action(action='remesh',
                                image=ImageResource('remesh.png'),
                                style='push',
                                tooltip = 'Incremental remesh',
                                show_label=False)

    move_vertices = Action(action='move_vertices',
                           image=ImageResource('movevertex.png'),
                           style='toggle',
                           tooltip = 'Move vertices',
                           show_label=False)

    tool_bar = ToolBar(flip_edges,
                       split_edges,
                       collapse_edges,
                       catmull_clark,
                       loop_subdivision,
                       incremental_remesh,
                       move_vertices,
                       show_labels=False)

    #--------------------------------------------------------------------------
    view = View(
                Item('scene',
                     editor=editor,
                     show_label=False,
                     resizable=True,
                     height=800,
                     width=1200,
                     ),
                resizable=True,
                title = 'Mesh Check',
                icon = ImageResource('logo.png'),
                toolbar = tool_bar
                )
    #--------------------------------------------------------------------------

    def __init__(self, mesh):
        HasTraits.__init__(self)
        self.scenemanager = SceneManager(scene_model=self.scene)
        self.scenemanager.add(mesh, rename=True)
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))
        self.mesh_manager.add_plot_callback(self.update_plot)

    @property
    def mesh(self):
        return self.mesh_manager.mesh

    @property
    def mesh_manager(self):
        return self.scenemanager.current_object

    def start(self):
        self.configure_traits()

    def update_plot(self):
        self.mesh_manager.plot_edges(color=(0.6,0.6,0.6))
        self.mesh_manager.plot_faces(color=(0.6,0.6,0.6))

    def _set_state(self, name):
        if name != 'flip_edges':
            self.flip_edges.checked = False
        if name != 'split_edges':
            self.split_edges.checked = False
        if name != 'collapse_edges':
            self.collapse_edges.checked = False
        if name != 'move_vertices':
            self.move_vertices.checked = False

    #--------------------------------------------------------------------------

    def flip(self, edge_index):
        self.mesh.flip_edge(edge_index)
        self.update_plot()

    def collapse(self, edge_index):
        self.mesh.collapse_edge(edge_index)
        self.mesh_manager.update_plot()

    def split(self, edge_index):
        self.mesh.split_edge(edge_index)
        self.mesh_manager.update_plot()

    def flip_edges(self):
        if self.flip_edges.checked:
            self._set_state('flip_edges')
            self.mesh_manager.on_edge_selection(self.flip)
        else:
            self.mesh_manager.select_edges_off()

    def collapse_edges(self):
        if self.collapse_edges.checked:
            self._set_state('collapse_edges')
            self.mesh_manager.on_edge_selection(self.collapse)
        else:
            self.mesh_manager.select_edges_off()

    def split_edges(self):
        if self.split_edges.checked:
            self._set_state('split_edges')
            self.mesh_manager.on_edge_selection(self.split)
        else:
            self.mesh_manager.select_edges_off()

    def catmull_clark(self):
        self._set_state(None)
        self.mesh.catmull_clark()
        self.mesh_manager.update_plot()

    def loop(self):
        self._set_state(None)
        self.mesh.loop()
        self.mesh_manager.update_plot()

    def remesh(self):
        self._set_state(None)
        self.remesher.mesh_manager = self.mesh_manager
        self.remesher.start()

    def move_vertices(self):
        if self.move_vertices.checked:
            self._set_state('move_vertices')
            def callback():
                self.mesh_manager.update_plot()
            self.mesh_manager.move_vertices(callback)
        else:
            self.mesh_manager.move_vertices_off()




# -----------------------------------------------------------------------------

def toolbar(mesh):
    m = MeshEditToolbar(mesh)
    m.start()
