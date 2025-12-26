#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from traits.api import HasTraits, Instance, Button, on_trait_change, Bool

from traitsui.api import View, Item, HSplit, VGroup

from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel

from mayavi.core.ui.mayavi_scene import MayaviScene

from pyface.image_resource import ImageResource

#------------------------------------------------------------------------------

from geometrylab.gui.scenemanager import SceneManager

# -----------------------------------------------------------------------------

'''check.py: an interactive checker for mesh connectivity and normals'''

__author__ = 'Davide Pellis'


###############################################################################
#################################### VIEWER ###################################
###############################################################################

class Selector(HasTraits):

    scene = Instance(MlabSceneModel, ())

    editor = SceneEditor(scene_class=MayaviScene)

    select_vertices = Bool()

    ok = Button()

    #--------------------------------------------------------------------------
    view = View(
                HSplit(VGroup(
                              Item('select_vertices'),
                              Item('ok'),
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
        if name != 'select_vertices':
            self.half_edge = False

    #--------------------------------------------------------------------------



    @on_trait_change('select_vertices')
    def select_vertices_check(self):
        if self.select_vertices:
            self.set_state('select_vertices')
            self.mesh_manager.select_vertices()
        else:
            self.mesh_manager.select_vertices_off()





# -----------------------------------------------------------------------------

def select_vertices(mesh):
    m = Selector(mesh)
    m.start()
    return m.mesh_manager.selected_vertices
