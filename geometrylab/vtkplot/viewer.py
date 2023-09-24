# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from traits.api import HasTraits, Instance, Property, Enum, Button,String,\
                       on_trait_change, Float, Bool, Int,Constant, ReadOnly,\
                       List, Array

from traitsui.api import View, Item, HSplit, VSplit, InstanceEditor, HGroup,\
                         Group, ListEditor, Tabbed, VGroup, CheckListEditor,\
                         ArrayEditor

from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel



from mayavi.core.ui.mayavi_scene import MayaviScene

from pyface.image_resource import ImageResource

import os

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.gui.geolabscene import GeolabScene

from geometrylab.vtkplot.plotmanager import PlotManager

# -----------------------------------------------------------------------------

'''viewer.py: The viewer for vtkplot souce classes'''

__author__ = 'Davide Pellis'

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                    Viewer
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class Viewer(HasTraits):

    # -------------------------------------------------------------------------
    #                                  Traits
    # -------------------------------------------------------------------------

    scene = Instance(MlabSceneModel, ())

    editor = SceneEditor(scene_class=GeolabScene)

    # -------------------------------------------------------------------------
    #                                   View
    # -------------------------------------------------------------------------

    view = View(Item('scene',
                     editor=editor,
                     show_label=False,
                     resizable=True,
                     height=800,
                     width=1200,
                     ),
                resizable=True,
                title = 'Geometrylab Viewer',
                icon = ImageResource(path + '/gui/img/new2/logo3.png')
                )

    # -------------------------------------------------------------------------
    #                                 Initialize
    # -------------------------------------------------------------------------

    def __init__(self, objects):
        HasTraits.__init__(self)

        self._plotmanager = PlotManager(scene_model=self.scene)

        self._generate_data(objects)

        self._position = None

    @property
    def position(self):
        return self.editor.scene_class._position

    @position.setter
    def position(self, position):
        self.editor.scene_class._position = position

    @property
    def background(self):
        return self._plotmanager.background

    @background.setter
    def background(self, background):
        self._plotmanager.background = background

    # -------------------------------------------------------------------------
    #                                 Methods
    # -------------------------------------------------------------------------

    def start(self):
        self.configure_traits()

    def _generate_data(self, objects):
        self._plotmanager.add(objects, rename=True)



# -----------------------------------------------------------------------------
#                                 View Function
# -----------------------------------------------------------------------------

def view(objects, position=None, background=(1,1,1)):
    viewer = Viewer(objects)
    viewer.background = background
    viewer.position = position
    viewer.start()
