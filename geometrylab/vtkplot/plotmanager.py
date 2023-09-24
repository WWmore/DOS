#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

from tvtk.api import tvtk

from mayavi.core.api import ModuleManager

from traits.api import HasTraits

from geometrylab.vtkplot.pointsource import Points

from geometrylab.vtkplot.facesource import Faces

from geometrylab.vtkplot.edgesource import Edges

from geometrylab.vtkplot.vectorsource import Vectors

from geometrylab.vtkplot.polylinesource import Polyline

#------------------------------------------------------------------------------

'''plotmanager.py: The plot manager class'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Plot Manager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class PlotManager(HasTraits):

    def __init__(self, scene_model=None, name='plot_0'):
        HasTraits.__init__(self)

        self.name = name

        self._scene_model = scene_model

        self._sources = {}

        self._widgets = {}

        self.picker_tolerance = 0.003

        self._position = None

        self._plot_callbacks = {}

        self._cross_callbacks = {}

    # -------------------------------------------------------------------------
    #                                Properties
    # -------------------------------------------------------------------------

    @property
    def type(self):
        return 'Plot_manager'

    @property
    def scene_model(self):
        return self._scene_model

    @scene_model.setter
    def scene_model(self, scene_model):
        self._scene_model = scene_model

    @property
    def engine(self):
        return self.scene_model.engine

    @property
    def scene(self):
        return self._scene_model.mayavi_scene

    @property
    def background(self):
        return self.scene_model.background

    @background.setter
    def background(self, color):
        self.scene_model.background = color

    @property
    def position(self):
        return self.get_position()

    @position.setter
    def position(self, position):
        if position is not None:
            #self.set_view(position)
            self._position = position

    @property
    def sources(self):
        return self._sources

    # -------------------------------------------------------------------------
    #                                    view
    # -------------------------------------------------------------------------

    def parallel_projection(self, parallel=True):
        self.scene.scene.camera.parallel_projection = parallel

    def z_view(self):
        self.scene.scene.z_plus_view()
        self.scene.scene.camera.compute_view_plane_normal()

    def interaction_2d(self):
        self.scene.scene.interactor.interactor_style = tvtk.InteractorStyleImage()

    # -------------------------------------------------------------------------
    #                            Pipeline management
    # -------------------------------------------------------------------------

    def add_widget(self, widget, name='widget'):
        self._widgets[name] = widget
        self.scene_model.add_actors(widget)

    def remove_widgets(self):
        for key in self._widgets:
            try:
                self.scene_model.remove_actor(self._widgets[key])
            except:
                pass
        self._widgets = {}

    def add(self, objects, plot=True, pickable=False, rename=False):
        if not isinstance(objects, list):
            objects = [objects]
        e = self.scene_model.engine
        self.fix_view()
        n = 0
        for obj in objects:
            if rename:
                obj.name = 'obj_' + str(n)
                n += 1
            try:
                self._sources[obj.name].source.remove()
            except:
                pass
            src = obj.source
            self._sources[obj.name] = obj
            if plot:
                e.add_source(src, scene=self.scene_model)
                obj.on_scene = True
            if not pickable:
                self.set_pickable(pick=False, name=obj.name)
        self.apply_view()

    def clear(self):
        self.fix_view()
        for key in self._sources:
            try:
                self._sources[key].source.remove()
                self._sources[key].on_scene = False
            except:
                pass
        self._sources = {}
        self.remove_widgets()
        self.apply_view()

    def hide(self, name=None):
        self.fix_view()
        if name is None:
            for key in self._sources:
                try:
                    self._sources[key].source.remove()
                    self._sources[key].on_scene = False
                except:
                    pass
        else:
            try:
               self._sources[name].source.remove()
               self._sources[name].on_scene = False
            except:
                pass
        self.apply_view()

    def remove(self, names):
        if type(names) is str:
            names = [names]
        self.fix_view()
        for key in names:
            try:
                self._sources[key].source.remove()
                del self._sources[key]
            except:
                pass
        self.apply_view()

    def keep_only(self, keep):
        for key in self._sources:
            if key not in keep:
                try:
                    self._sources[key].remove()
                    self.remove(key)
                except:
                    pass

    def show(self, name):
        e = self.scene_model.engine
        self.fix_view()
        try:
            e.add_source(self._sources[name].source, scene=self.scene_model)
            self._sources[name].on_scene = True
        except:
            pass
        self.apply_view()

    def update(self, name, **kwargs):
        self.fix_view()
        self._sources[name].update(**kwargs)
        if not self._sources[name].on_scene:
            e = self.scene_model.engine
            e.add_source(self._sources[name].source, scene=self.scene_model)
            self._sources[name].on_scene = True
        self.apply_view()

    def has_source(self, name):
        if name in self._sources:
            return True
        else:
            return False

    def is_on_scene(self, name):
        if name in self._sources:
            if self._sources[name].on_scene:
                return True
        return False

    def get_source(self, name):
        try:
            return self.sources[name]
        except:
            return None

    # -------------------------------------------------------------------------
    #                                 Output
    # -------------------------------------------------------------------------

    def save(self, name, size=None):
        scene = self.scene_model.engine.current_scene.scene
        #scene.render_window.point_smoothing = True
        #scene.render_window.line_smoothing = True
        #scene.render_window.polygon_smoothing = True
        #scene.render_window.multi_samples = 8
        #scene.render_window.aa_frames = 16
        scene.save_png(name)
        #self.scene_model.save(name, size)

    # -------------------------------------------------------------------------
    #                              View management
    # -------------------------------------------------------------------------

    def disable_render(self):
        self.scene.scene.disable_render = True

    def enable_render(self):
        self.scene.scene.disable_render = False

    def fix_view(self):
        try:
            self.scene.scene.disable_render = True
            cc = self.scene.scene.camera
            self.__orig_pos = cc.position
            self.__orig_fp = cc.focal_point
            self.__orig_view_angle = cc.view_angle
            self.__orig_view_up = cc.view_up
            self.__orig_clipping_range = cc.clipping_range
            self.__parallel_projection = cc.parallel_projection
        except:
            pass

    def apply_view(self):
        try:
            cc = self.scene.scene.camera
            cc.position = self.__orig_pos
            cc.focal_point = self.__orig_fp
            cc.view_angle =  self.__orig_view_angle
            cc.view_up = self.__orig_view_up
            cc.clipping_range = self.__orig_clipping_range
            cc.parallel_projection = self.__parallel_projection
            self.scene.scene.disable_render=False
        except:
            pass

    def get_position(self):
        scene = self.scene.scene
        cc = scene.camera
        p = [cc.position[0],cc.position[1],cc.position[2]]
        p.extend([cc.view_up[0],cc.view_up[1],cc.view_up[2]])
        p.extend([cc.focal_point[0],cc.focal_point[1],cc.focal_point[2]])
        p.extend([cc.view_angle])
        p.extend([cc.clipping_range[0],cc.clipping_range[1]])
        p.extend([scene.parallel_projection])
        return p

    def camera_zRotation(self, center, degrees):
        position = np.array(self.get_position())
        P = position[0:2]
        F = position[6:8]
        cos = np.cos(np.radians(degrees))
        sin = np.sin(np.radians(degrees))
        V = P - center[0:2]
        Vx = V[0]*cos - V[1]*sin + center[0]
        Vy = V[0]*sin + V[1]*cos + center[1]
        position[0] = Vx
        position[1] = Vy
        V = F - center[0:2]
        Vx = V[0]*cos - V[1]*sin + center[0]
        Vy = V[0]*sin + V[1]*cos + center[1]
        position[6] = Vx
        position[7] = Vy
        position[3:6] = np.array([0,0,1])
        self.set_view(position)

    def save_view(self, print_position=True):
        self._position = self.get_position()
        if print_position:
            print('position = ' + str(self._position))

    def set_view(self, position=None):
        cc = self.scene.scene.camera
        p = self._position
        if position is not None:
            p = position
        if p is not None:
            cc.position = np.array([p[0],p[1],p[2]])
            cc.view_up = np.array([p[3],p[4],p[5]])
            cc.focal_point = np.array([p[6],p[7],p[8]])
            cc.view_angle = np.float(p[9])
            cc.clipping_range = np.array([p[10],p[11]])
            try:
                cc.parallel_projection = p[12]
            except:
                pass
            self.scene.render()

    def render(self):
        self.scene.scene._renwin.render()

    # -------------------------------------------------------------------------
    #                            Get from pipeline
    # -------------------------------------------------------------------------

    def get_actor(self, name):
        stop = False
        try:
            obj = self._sources[name].source
        except KeyError:
            return
        while not stop:
            obj = obj.children[0]
            if isinstance(obj, ModuleManager):
                stop = True
        actor = obj.children[0].actor
        return actor

    # -------------------------------------------------------------------------
    #                                 Picker
    # -------------------------------------------------------------------------

    def set_pickable(self, pick=True, name=None):
        from mayavi.core.api import ModuleManager
        if name is None:
            for key in self._sources:
                obj = self._sources[key].source
                if pick:
                    pick = 1
                else:
                    pick = 0
                stop = False
                while not stop:
                    obj = obj.children[0]
                    if isinstance(obj, ModuleManager):
                        stop = True
                actor = obj.children[0].actor.actors[0]
                actor.pickable = pick
        else:
            try:
                obj = self._sources[name].source
            except:
                return
            if pick:
                pick = 1
            else:
                pick = 0
            stop = False
            while not stop:
                obj = obj.children[0]
                if isinstance(obj, ModuleManager):
                    stop = True
            actor = obj.children[0].actor.actors[0]
            actor.pickable = pick

    def picker_callback(self, routine, mode='point', name=None, add=False):
        if name is not None:
            self.set_pickable(False)
            self.set_pickable(True, name)
        else:
            self.set_pickable(True)
        def picker_callback(picker):
            pick_id = picker.cell_id
            if mode == 'point':
                pick_id = picker.cell_id
            routine(pick_id)
        #s = self.scene_model.engine.current_scene
        s = self.scene
        p = s._mouse_pick_dispatcher.callbacks
        if add:
            s.on_mouse_pick(picker_callback,'cell')
        elif len(p) == 0:
            s.on_mouse_pick(picker_callback,'cell')
        else:
            p[0] = (picker_callback,'cell','Left')
        a = s._mouse_pick_dispatcher._active_pickers['cell']
        a.tolerance = self.picker_tolerance

    def picker_off(self):
        #scene = self.scene_model.engine.current_scene
        p = self.scene._mouse_pick_dispatcher.callbacks
        def picker_callback(picker):
            return
        for i in range(len(p)):
            p[i] = (picker_callback,'cell','Left')

    def iterate(self, function, times=1):
        for i in range(times):
            function()
            #self.scene_model.engine.current_scene.scene._renwin.render()
            self.scene.scene._renwin.render()

    # -------------------------------------------------------------------------
    #                                 Plot
    # -------------------------------------------------------------------------

    def update_plot(self, main=True):
        self.disable_render()
        for key in self._plot_callbacks:
            callback = self._plot_callbacks[key]
            callback()
        if main:
            for key in self._cross_callbacks:
                callback = self._cross_callbacks[key]
                callback(main=False)
        self.enable_render()

    def add_plot_callback(self, callback, name=None):
        if callable(callback):
            if name is None:
                name = ('plot_callback_{}').format(self.__counter)
                self.__counter += 1
            self._plot_callbacks[name] = callback

    def add_cross_callback(self, callback, name=None, **kwargs):
        if callable(callback):
            if name is None:
                name = ('cross_callback_{}').format(self.__counter)
                self.__counter += 1
            self._cross_callbacks[name] = callback

    def remove_plot_callback(self, name):
        try:
            del(self._plot_callbacks[name])
        except:
            pass

    # -------------------------------------------------------------------------
    #                              Plot Functions
    # -------------------------------------------------------------------------

    def plot_points(self, points, **kwargs):
        name = kwargs.pop('name', 'points')
        if name in self._sources:
            kwargs['points'] = points
            self.update(name, **kwargs)
        else:
            kwargs['name'] = name
            P = Points(points, **kwargs)
            self.add(P)

    def plot_polyline(self, polyline, **kwargs):
        name = kwargs.pop('name', 'polyline')
        if name in self._sources:
            kwargs['polyline'] = polyline
            self.update(name, **kwargs)
        else:
            kwargs['name'] = name
            P = Polyline(polyline, **kwargs)
            self.add(P)

    def plot_mesh_faces(self, mesh, **kwargs):
        name = kwargs.pop('name', 'mesh-faces')
        if name in self._sources:
            kwargs['mesh'] = mesh
            self.update(name, **kwargs)
        else:
            kwargs['name'] = name
            P = Faces(mesh, **kwargs)
            self.add(P)

    def plot_mesh_edges(self, mesh, **kwargs):
        name = kwargs.pop('name', 'mesh-edges')
        if name in self._sources:
            kwargs['mesh'] = mesh
            self.update(name, **kwargs)
        else:
            kwargs['name'] = name
            P = Edges(mesh, **kwargs)
            self.add(P)

    def plot_vectors(self, vectors, **kwargs):
        name = kwargs.pop('name', 'vectors')
        if name in self._sources:
            kwargs['vectors'] = vectors
            self.update(name, **kwargs)
        else:
            kwargs['name'] = name
            P = Vectors(vectors, **kwargs)
            self.add(P)

