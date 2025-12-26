#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

# -----------------------------------------------------------------------------

from geometrylab.gui.scenemanager import SceneManager

# -----------------------------------------------------------------------------

'''plotmanager.py: The plot manager class'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                            Multiple Scene Manager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class MultiSceneManager(SceneManager):

    def __init__(self, scene_models=None):
        SceneManager.__init__(self)

        self._scene_managers = {}

        self._object_scenes = {}

        self._scene_register = {}

        #self._engine = Engine()

        #self.engine.start()

    @property
    def object_keys(self):
        return list(self._object_scenes.keys())

    @property
    def scene_keys(self):
        return list(self._scene_managers.keys())

    @property
    def engine(self):
        #scenemodel = self._scene_managers[self.scene_keys[0]].scene_model
        #return scenemodel.engine
        return None

    @property
    def scene_managers(self):
        return self._scene_managers

    @property
    def current_object(self):
        if self._selected is None:
            try:
                return self.get_object(self.object_keys[0])
            except:
                return None
        else:
            return self.get_object(self._selected)

    @property
    def last_object(self):
        if self._last_object is None:
            try:
                return self.get_object(self.object_keys[0])
            except:
                return None
        else:
            return self.get_object(self._last_object)

    @property
    def current_object_type(self):
        if self._selected is None:
            try:
                return self.get_object(self.object_keys[0]).type
            except:
                return 'none'
        else:
            return self.get_object(self._selected).type

    #--------------------------------------------------------------------------
    #                                 Scenes
    #--------------------------------------------------------------------------

    def clear_scenes(self):
        for scene in self.engine.scenes:
            scene.remove()

    def add_scene_model(self, scene_model, name):
        manager = SceneManager(scene_model)
        self._scene_managers[name] = manager
        self._scene_register[scene_model.mayavi_scene] = manager

    def get_scene(self, name):
        return self._scene_managers[name]

    def _format_scene(self, scene):
        if type(scene) is dict:
            scene = scene.get('scene', None)
        if scene is None:
            scene = 'scene_0'
        return scene

    def set_background(self, color):
        for key in self._scene_managers:
            self._scene_managers[key].background = color

    #--------------------------------------------------------------------------
    #                                 Objects
    #--------------------------------------------------------------------------

    def add_object(self, geometry, name=None, scene=None):
        if name is None:
            name = ('obj_{}').format(self._counter)
            self._counter += 1
        scene = self._format_scene(scene)
        self.get_scene(scene).add_object(geometry, name=name)
        self._last_object = name
        obj = self.get_scene(scene).get_object(name)
        self._object_scenes[name] = scene
        self._objects[name] = obj
        for key in self._objects:
            other_obj = self._objects[key]
            if obj != other_obj:
                if (id(other_obj.geometry.vertices) ==
                                        id(obj.geometry.vertices)):
                    obj.add_cross_callback(other_obj.update_plot, name=name)
                    other_obj.add_cross_callback(obj.update_plot, name=name)

    def add(self, objects, scene=None, rename=False):
        scene = self._format_scene(scene)
        self.get_scene(scene).add(objects, rename=rename)
        for scene in self.scene_managers:
            for key, obj in self.scene_managers[scene]._objects.items():
                self._object_scenes[key] = scene
                self._objects[key] = obj

    def get_object(self, name):
        scene = self.get_scene(self._object_scenes[name])
        obj = scene.get_object(name)
        return obj

    def clear(self, scene=None):
        scene = self._format_scene(scene)
        self.get_scene(scene).clear()

    def hide(self, name):
        for key in self.scene_managers:
            self.get_scene(key).hide(name)

    #--------------------------------------------------------------------------
    #                                 Plot
    #--------------------------------------------------------------------------

    def disable_render(self):
        for key in self.scene_managers:
            self.get_scene(key).disable_render()

    def enable_render(self):
        for key in self.scene_managers:
            self.get_scene(key).enable_render()

    def plot_points(self, points, **kwargs):
        scene = self._format_scene(kwargs)
        self.get_scene(scene).plot_points(points, **kwargs)

    def plot_polyline(self, polyline, **kwargs):
        scene = self._format_scene(kwargs)
        self.get_scene(scene).plot_polyline(polyline, **kwargs)

    def update_plot(self, scene=None):
        self.disable_render()
        if scene is None:
            for key in self._scene_managers:
                self.get_scene(key).update_plot()
        else:
            scene = self._format_scene(scene)
            self.get_scene(scene).update_plot()
        for key in self._plot_callbacks:
            self._plot_callbacks[key]()
        self.enable_render()

    #--------------------------------------------------------------------------
    #                                 Select
    #--------------------------------------------------------------------------

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
        for key in self.scene_managers:
            s = self.scene_managers[key].scene_model.mayavi_scene
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
            obj.highlight_off()
        def picker_callback(picker):
            return
        for key in self.scene_managers:
            s = self.scene_managers[key].scene_model.mayavi_scene
            p = s._mouse_pick_dispatcher.callbacks
            for i in range(len(p)):
                p[i] = (picker_callback,'cell','Left')
