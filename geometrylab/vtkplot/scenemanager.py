#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

# -----------------------------------------------------------------------------

from geometrylab.vtkplot.plotmanager import PlotManager

from geometrylab.vtkplot.meshplotmanager import MeshPlotManager


# -----------------------------------------------------------------------------

'''plotmanager.py: The plot manager class'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Plot Manager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class SceneManager(PlotManager):

    def __init__(self, scene_model=None):
        super(SceneManager, self).__init__(scene_model)

        self._objects = {}

        self._selected = None

        self._counter = 0

    # -------------------------------------------------------------------------
    #                                Properties
    # -------------------------------------------------------------------------

    @property
    def scene_model(self):
        return self._scene_model

    @scene_model.setter
    def scene_model(self, scene_model):
        self._scene_model = scene_model
        for obj in self._objects:
            obj.scene_model = scene_model

    @property
    def current_object(self):
        if self._selected is None:
            try:
                return self._objects[self._objects.keys()[0]]
            except:
                return None
        return self._objects[self._selected]

    # -------------------------------------------------------------------------
    #                            Get from pipeline
    # -------------------------------------------------------------------------

    def add(self, objects, plot=True, rename=False, name=None):
        if not isinstance(objects, list):
            objects = [objects]
        for obj in objects:
            if obj.type is 'Mesh':
                pm = MeshPlotManager(self.scene_model)
                pm.mesh = obj
                obj = pm
            if name is not None:
                obj.name = name
            elif rename:
                obj.name = ('obj_{}').format(self._counter)
                self._counter += 1
            if obj.type is 'Mesh_plot_manager':
                self._objects[obj.name] = obj
                obj.scene_model = self.scene_model
                if plot:
                    obj.update_plot()
            else:
                super(SceneManager, self).add(obj)

    def remove(self, names):
        if type(names) is str:
            names = [names]
        for key in names:
            try:
                self._objects[key].clear()
                del self._objects[key]
            except:
                pass
        super(SceneManager, self).remove(names)

    def update_plot(self, **kwargs):
        for key in self._objects:
            obj = self._objects[key]
            obj.update_plot()
        super(SceneManager, self).update(*kwargs)

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

    # -------------------------------------------------------------------------
    #                                  Selection
    # -------------------------------------------------------------------------

    def select_object(self):
        self.current_object.hiligth()
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
                point_id = picker.point_id
                pp = picker.mapper_position
                d_min = float('inf')
                selected = None
                for key in self._objects:
                    try:
                        obj = self._objects[key]
                        p = obj.mesh.vertices[point_id]
                        d = (pp[0]-p[0])**2 + (pp[1]-p[1])**2 + (pp[2]-p[2])**2
                        if d < d_min:
                            d_min = d
                            selected = key
                    except IndexError:
                        pass
                self.current_object.hiligth_off()
                self._selected = selected
                self.current_object.hiligth()
        s = self.scene_model.engine.current_scene
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
            obj.selection_off()
            obj.hiligth_off()
        scene = self.scene_model.engine.current_scene
        p = scene._mouse_pick_dispatcher.callbacks
        def picker_callback(picker):
            return
        for i in range(len(p)):
            p[i] = (picker_callback,'cell','Left')