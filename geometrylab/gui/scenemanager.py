#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

# -----------------------------------------------------------------------------

from geometrylab.vtkplot.plotmanager import PlotManager

from geometrylab.gui.geolabmesh import GeolabMesh

from geometrylab.gui.geolabpoints import GeolabPoints

from geometrylab.optimization.gridshell import Gridshell

# -----------------------------------------------------------------------------

'''plotmanager.py: The plot manager class'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                Scene Manager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class SceneManager(PlotManager):

    def __init__(self, scene_model=None):
        PlotManager.__init__(self, scene_model)

        self._objects = {}

        self._selected = None

        self._last_object = None

        self._counter = 0

        self._handler = None

        self._geometry_change_callbacks = []

        self._object_changed_callbacks = []

        self._object_change_callbacks = []

        self._object_save_callbacks = []

        self._object_open_callbacks = []

    # -------------------------------------------------------------------------
    #                                Properties
    # -------------------------------------------------------------------------

    @property
    def scene_model(self):
        return self._scene_model

    @scene_model.setter
    def scene_model(self, scene_model):
        self._scene_model = scene_model
        self._scene = scene_model.mayavi_scene
        for obj in self._objects:
            obj.scene_model = scene_model

    @property
    def current_object(self):
        if self._selected is None:
            try:
                return self._objects[self._objects.keys()[0]]
            except:
                return None
        try:
            return self._objects[self._selected]
        except:
            self._selected = None
            return self.current_object

    @property
    def last_object(self):
        if self._last_object is None:
            try:
                return self._objects[self._objects.keys()[0]]
            except:
                return None
        try:
            return self._objects[self._last_object]
        except:
            self._last_object = None
            return self.last_object

    def add_callback(self, key, callback):
        if key == 'object_change':
            self._object_change_callbacks.append(callback)
        elif key == 'object_changed':
            self._object_changed_callbacks.append(callback)
        elif key == 'object_save':
            self._object_save_callbacks.append(callback)
        elif key == 'object_open':
            self._object_open_callbacks.append(callback)

    def object_changed(self):
        for callback in self._object_changed_callbacks:
            callback()

    def object_change(self):
        for callback in self._object_change_callbacks:
            callback()

    def object_save(self, file_name):
        for callback in self._object_save_callbacks:
            callback(file_name)

    def object_open(self, file_name, geometry):
        for callback in self._object_open_callbacks:
            callback(file_name, geometry)

    # -------------------------------------------------------------------------
    #                            Get from pipeline
    # -------------------------------------------------------------------------

    def clear_scene(self):
        for scene in self.engine.scenes:
           if scene != self.scene_model:
               pass#scene.remove()

    def add_object(self, geometry, name=None):
        if geometry.type == 'Mesh':
            if not hasattr(geometry, 'constrained_vertices'):
                g = Gridshell()
                g.import_mesh(geometry)
                geometry = g
            GM = GeolabMesh(self.scene_model)
            GM.geometry = geometry
            obj = GM
        if geometry.type == 'Points':
            GP = GeolabPoints(self.scene_model)
            GP.geometry = geometry
            obj = GP
        if name is not None:
            obj.name = name
        else:
            obj.name = ('obj_{}').format(self._counter)
            self._counter += 1
        if len(self._objects) == 0:
            self._selected = obj.name
        self._last_object = obj.name
        self.remove(obj.name)
        self._objects[obj.name] = obj
        obj.scene_model = self.scene_model

    def remove(self, names=None):
        if names is None:
            if self._selected is None:
                names = [self._objects.keys()[0]]
            else:
                names = [str(self._selected)]
        if type(names) is str:
            names = [names]
        for key in names:
            try:
                self._objects[key].clear()
                del self._objects[key]
            except:
                pass
        PlotManager.remove(self, names)

    def update_plot(self):
        PlotManager.update_plot(self)
        for key in self._objects:
            obj = self._objects[key]
            obj.update_plot()

    def hide(self, name=None):
        if name is None:
            super(SceneManager, self).hide()
            for key in self._objects:
                try:
                    self._objects[key].hide()
                except:
                    pass
        else:
            super(SceneManager, self).hide(name)
            try:
                self._objects[name].hide()
            except:
                pass

    def get_object(self, name):
        try:
            obj = self._objects[name]
        except:
            obj = None
        return obj

    # -------------------------------------------------------------------------
    #                                  Selection
    # -------------------------------------------------------------------------

    def select_object(self):
        self.current_object.highlight()
        ranges = {}
        count = 0
        for key in self._objects:
            obj = self._objects[key]
            N = obj.selection_on()
            ranges[key] = [count, count+N]
            count += N
        def picker_callback(picker):
            pick_id = picker.cell_id
            if pick_id >= 0:
                for key in self._objects:
                    try:
                        obj = self._objects[key]
                        if picker.actor in obj.object_selection_actors:
                            self.current_object.highlight_off()
                            self.object_change()
                            self._selected = key
                            self.object_changed()
                            self.current_object.highlight()
                    except IndexError:
                        pass

        s = self.engine.current_scene
        p = s._mouse_pick_dispatcher.callbacks
        if len(p) == 0:
            s.on_mouse_pick(picker_callback,'cell')
        else:
            p[0] = (picker_callback,'cell','Left')
        a = s._mouse_pick_dispatcher._active_pickers['cell']
        a.tolerance = self.picker_tolerance

    def select_off(self):
        for key in self._objects:
            obj = self._objects[key]
            obj.select_off()
            obj.hilight_off()
        scene = self.engine.current_scene
        p = scene._mouse_pick_dispatcher.callbacks
        def picker_callback(picker):
            return
        for i in range(len(p)):
            p[i] = (picker_callback,'cell','Left')