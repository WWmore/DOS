#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from traits.api import HasTraits, Button, String, on_trait_change, Float, Bool, Array, Range

from traitsui.api import View, Item,  HGroup, Group, VGroup,ArrayEditor, Controller

from pyface.image_resource import ImageResource

# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                Tools Handler
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class T_Handler(Controller):

    def close(self, info, is_ok):
        info.object._closed = True
        Controller.close(self ,info, is_ok)
        return True


class Tool(HasTraits):

    def __init__(self):
        HasTraits.__init__(self)
        self._scenemanager = None

    @property
    def mesh(self):
        return self.meshmanager.mesh

    @property
    def meshmanager(self):
        if self._scenemanager.current_object.type is not 'Mesh_plot_manager':
            return None
        return self._scenemanager.current_object

    @property
    def scenemanager(self):
        return self._scenemanager

    @property
    def geolab(self):
        return self._scenemanager

    @scenemanager.setter
    def scenemanager(self, scene_manager):
        self._scenemanager = scene_manager



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                              Incremental Remesh
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class Loads(Tool):

    _closed = Bool(True)

    _handler = T_Handler()

    area_load = Float(0,label='Area load')

    beam_load = Float(1,label='Beam load')

    vector = Array(np.float32,(1,3), editor = ArrayEditor(width=20))

    apply_force = Button(label='Apply force')

    reset_forces = Button(label='Reset forces')

    view = View(VGroup('area_load',
                       'beam_load',
                       'vector',
                       HGroup(
                              Item('apply_force',
                                   show_label=False,
                                   resizable=True),
                              Item('reset_forces',
                                   show_label=False,
                                   resizable=True),
                              ),
                      show_border=True),
                title='Loads',
                handler=_handler,
                width=300,
                icon=ImageResource('img/new/applyforce.png'))


    def start(self):
        if not self._closed:
            self._handler.info.ui.dispose()
        self.configure_traits()
        self._closed = False
        self.meshmanager.hide('forces')
        self.plot_loads()

    def close(self):
        try:
            self._closed = True
            self._handler.info.ui.dispose()
        except:
            pass

    @on_trait_change('_closed')
    def hide_loads(self):
        self.meshmanager.hide('apply_loads')
        self.meshmanager.update_plot()

    @on_trait_change('beam_load, area_load')
    def update_load(self):
        self.mesh.area_load = self.area_load
        self.mesh.beam_load = self.beam_load
        self.plot_loads()

    @on_trait_change('apply_force')
    def apply_force_fired(self):
        selected = self.meshmanager.selected_vertices
        self.mesh.apply_force(self.vector[0,:], selected)
        self.plot_loads()

    @on_trait_change('reset_forces')
    def reset_forces_fired(self):
        self.mesh.reset_forces()
        self.plot_loads()

    def plot_loads(self):
        self.meshmanager.plot_vectors(vectors = self.mesh.loads(),
                                     color = 'blue-red',
                                     glyph_type = '3D-arrow',
                                     position = 'head',
                                     name = 'apply_loads')


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                              Incremental Remesh
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class IncrementalRemesh(Tool):

    _closed = Bool(True)

    _handler = T_Handler()

    remesh_factor = Range(0.1, 1.5, 1)

    corner_tolerance = Range(-1.0, 1.0, 0.3)

    remesh_button = Button(label='Remesh')

    view = View(VGroup('remesh_factor',
                       'corner_tolerance',
                       Item('remesh_button',
                            show_label=False),
                       show_border=True),
                title='Incremental Remesh',
                handler=_handler,
                icon=ImageResource('img/new2/remesh.png'))

    def start(self):
        if not self._closed:
            self._handler.info.ui.dispose()
        self.configure_traits()
        self._closed = False

    def close(self):
        try:
            self._closed = True
            self._handler.info.ui.dispose()
        except:
            pass

    @on_trait_change('corner_tolerance')
    def set_corner_tolerance(self):
        self.mesh.corner_tolerance = self.corner_tolerance

    @on_trait_change('remesh_button')
    def incremental_remesh(self):
        if self.meshmanager is not None:
            self.mesh.incremental_remesh(self.remesh_factor)
            self.geolab.fix_edges_bug()
            self.meshmanager.update_plot()


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                              Corner Tolerance
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class CornerTolerance(Tool):

    corner_tolerance = Range(-1.0, 1.0, 0.3)

    _closed = Bool(True)

    _handler = T_Handler()

    view = View(Group('corner_tolerance',
                      show_border=True),
                title='Set corner tolerance',
                handler=_handler,
                icon=ImageResource('img/new2/corners.png'))


    @property
    def scenemanager(self):
        return self._scenemanager

    @scenemanager.setter
    def scenemanager(self, scene_manager):
        self._scenemanager = scene_manager
        ct = scene_manager.current_object.mesh.corner_tolerance
        self.corner_tolerance = ct

    def start(self):
        if not self._closed:
            self._handler.info.ui.dispose()
        self.configure_traits()
        self._closed = False
        self.plot_corners()
        self.scenemanager.current_object.add_plot_callback(self.plot_corners,
                                                        name='corners')

    def close(self):
        try:
            self._closed = True
            self._handler.info.ui.dispose()
        except:
            pass

    @on_trait_change('corner_tolerance')
    def set_corner_tolerance(self):
        self.mesh.corner_tolerance = self.corner_tolerance
        self.plot_corners()

    @on_trait_change('_closed')
    def hide_corners(self):
        self.scenemanager.current_object.hide('corners')
        self.scenemanager.current_object.remove_plot_callback('corners')

    def plot_corners(self):
        corners = self.mesh.mesh_corners()
        r = self.scenemanager.current_object.r
        self.scenemanager.current_object.plot_glyph(vertex_indices=corners,
                                                  glyph_type = 'sphere',
                                                  radius=2.2*r,
                                                  color='g',
                                                  name='corners')


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Save Mesh
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class SaveMesh(Tool):

    _closed = Bool(True)

    _handler = T_Handler()

    label = String('opt')

    save_button = Button(label='Save')

    view = View(HGroup('label',
                       Item('save_button',
                            show_label=False),
                       show_border=True),
                title='Save mesh',
                handler=_handler,
                icon=ImageResource('img/new/savemesh.png'))

    def start(self):
        if not self._closed:
            self._handler.info.ui.dispose()
        self.configure_traits()
        self._closed = False

    def close(self):
        try:
            self._closed = True
            self._handler.info.ui.dispose()
        except:
            pass

    @on_trait_change('save_button')
    def save_file(self):
        name = ('{}_{}').format(self.mesh.name, self.label)
        path = self.mesh.make_obj_file(name)
        self.scenemanager.object_save(path)
        if False:
            name = ('{}_stress_{}_{}').format(self.mesh.name, self.label,
                                            self.__counter)
            V1, V2 = self.mesh.equivalent_stress_directions()
            self.mesh.make_crossfield_obj_file(V1, V2, name)
        self.close()