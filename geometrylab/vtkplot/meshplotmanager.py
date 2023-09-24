#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

import numpy as np

#------------------------------------------------------------------------------

from geometrylab.vtkplot.plotmanager import PlotManager

from geometrylab.vtkplot.edgesource import Edges

from geometrylab.vtkplot.pointsource import Points

from geometrylab.vtkplot.facesource import Faces

from geometrylab.vtkplot.vectorsource import Vectors

#------------------------------------------------------------------------------

'''-'''

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Mesh PlotManager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class MeshPlotManager(PlotManager):

    def __init__(self, scene_model=None):
        PlotManager.__init__(self, scene_model)

        self.scale = 1

        self.glyph_scale = 1

        self.vector_scale = 1

        self._vertex_data = None

        self._face_data = None

        self._edge_data = None

        self.selected_vertices = []

        self.selected_faces = []

        self.selected_edges = []

        self.vertex_color = 'cornflower'

        self.edge_color = 'cornflower'

        self.face_color = 'cornflower'

        self.selection_color = 'yellow'

        self.view_mode = 'solid'

        self._show_virtual = False

        self._r = None

        self._g = None

        self._mesh = None

        self._face_sources = []

        self._edge_sources = []

        self._updating = True

        self._record = False

        self.__face_data = None

        self.__edge_data = None

        self.__dimensions = (0,0,0)

        self.__counter = 0

    # -------------------------------------------------------------------------
    #                              Properties
    # -------------------------------------------------------------------------

    @property
    def type(self):
        return 'Mesh_plot_manager'

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        if mesh.type is not 'Mesh':
            raise ValueError('*mesh* attribute must be a Mesh!')
        self._mesh = mesh
        self._set_r()
        self._set_g()
        self.vertex_data = np.zeros(mesh.V)
        self.face_data = np.zeros(mesh.F)
        self.edge_data = np.zeros(mesh.E)

    @property
    def geometry(self):
        return self._mesh

    @geometry.setter
    def geometry(self, mesh):
        if mesh.type is not 'Mesh':
            raise ValueError('*geometry* attribute must be a Mesh!')
        self.mesh = mesh

    @property
    def r(self):
        return self._r * self.scale

    @property
    def g(self):
        return self._r* self._g * self.glyph_scale

    @property
    def v(self):
        return 12 * self._r * self.vector_scale

    @property
    def updating(self):
        return self._updating

    @updating.setter
    def updating(self, updating):
        if type(updating) == bool:
            self._updating = updating
        else:
            raise ValueError('the updating attribute must be a boolean!')

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, bool):
        self._record = bool

    @property
    def vertex_data(self):
        if self._vertex_data is None:
            self._vertex_data = np.zeros(self.mesh.V)
        return self._vertex_data

    @vertex_data.setter
    def vertex_data(self, vertex_data):
        self._vertex_data = vertex_data

    @property
    def face_data(self):
        if self._face_data is None:
            self._face_data = np.zeros(self.mesh.F)
        return self._face_data

    @face_data.setter
    def face_data(self, face_data):
        self._face_data = face_data

    @property
    def edge_data(self):
        if self._edge_data is None:
            self._edge_data = np.zeros(self.mesh.E)
        return self._edge_data

    @edge_data.setter
    def edge_data(self, edge_data):
        self._edge_data = edge_data

    @property
    def object_selection_actors(self):
        return self.get_actor('virtual-mesh').actors

    # -------------------------------------------------------------------------
    #                              Size set
    # -------------------------------------------------------------------------

    def check_updates(self):
        updated = False
        if self.__dimensions[0] != self.mesh.V:
            updated = True
        if self.__dimensions[1] != self.mesh.E:
            updated = True
        if self.__dimensions[2] != self.mesh.F:
            updated = True
        if updated:
            self._set_r()
            self._set_g()
            self._vertex_data = None
            self._edge_data = None
            self._face_data = None
        self.__dimensions = (self.mesh.V, self.mesh.E, self.mesh.F)
        return updated

    def _set_r(self):
        if self._mesh is None or self._mesh.V == 0:
            self._r = None
        else:
            d = self._mesh.mean_edge_length()
            r = d/50 * max(1, np.log10(self._mesh.E))
            self._r = r

    def _set_g(self):
        if self._mesh is None or self._mesh.V == 0:
            self._g = None
        else:
            g = max(1, np.log10(self._mesh.E/100))
            self._g = g

    # -------------------------------------------------------------------------
    #                                Clear
    # -------------------------------------------------------------------------

    def clear(self, delete=True):
        super(MeshPlotManager, self).clear()
        self.picker_off()
        #self.remove_widgets()
        self.selected_vertices = []
        self.selected_edges = []
        self.selected_faces = []


    # -------------------------------------------------------------------------
    #                            Plot functions
    # -------------------------------------------------------------------------

    def initialize_plot(self):
        pass
        #self.plot_edges()

    def update_plot(self, **kwargs):
        self.check_updates()
        PlotManager.update_plot(self, **kwargs)
        self.plot_selected_vertices()

    def plot_faces(self, **kwargs):
        clear = False
        if 'name' not in kwargs:
            kwargs['name'] = 'faces'
        else:
            self._face_sources.append(kwargs['name'])
        if 'color' not in kwargs:
            kwargs['color'] = self.face_color
        if 'face_data' in kwargs:
            if self.__face_data is not 'face':
                clear = True
            if self.record:
                kwargs['face_data'] = self.face_data
            else:
                self.face_data = kwargs['face_data']
            self.__face_data = 'face'
        elif 'vertex_data' in kwargs:
            if self.__face_data is not 'vertex':
                clear = True
            if self.record:
                kwargs['vertex_data'] = self.vertex_data
            else:
                self.veretx_data = kwargs['vertex_data']
            self.__face_data = 'vertex'
        if kwargs['name'] not in self.sources or not self.updating or clear:
            F = Faces(self.mesh, **kwargs)
            self.add(F)
        else:
            self.update(**kwargs)

    def plot_edges(self, **kwargs):
        clear = False
        if 'update' not in kwargs:
            kwargs['update'] = True
        if 'name' not in kwargs:
            kwargs['name'] = 'edges'
        else:
            self._edge_sources.append(kwargs['name'])
        if 'tube_radius' not in kwargs:
            if self.view_mode is 'wireframe':
                tube_radius = None
            elif self.view_mode == '3d':
                tube_radius = None
            else:
                tube_radius = self.r
            kwargs['tube_radius'] = tube_radius
        if 'color' not in kwargs:
            kwargs['color'] = self.edge_color
        if 'edge_data' in kwargs:
            if self.__edge_data is not 'edge':
                clear = True
            if self.record:
                kwargs['edge_data'] = self.edge_data
            else:
                self.edge_data = kwargs['edge_data']
            self.__edge_data = 'edge'
        elif 'vertex_data' in kwargs:
            if self.__edge_data is not 'vertex':
                clear = True
            if self.record:
                kwargs['vertex_data'] = self.vertex_data
            else:
                self.veretx_data = kwargs['vertex_data']
            self.__edge_data = 'vertex'
        if kwargs['name'] not in self.sources or not self.updating or clear:
            M = Edges(self.mesh, **kwargs)
            self.add(M)
        else:
            self.update(**kwargs)

    def plot_vertices(self, **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = 'vertices'
        if 'radius' not in kwargs:
            kwargs['radius'] = self.r * 2.5
        if 'color' not in kwargs:
            kwargs['color'] = self.vertex_color
        if kwargs['name'] not in self.sources or not self.updating:
            P = Points(self.mesh, **kwargs)
            self.add(P)
        else:
            self.update(**kwargs)

    def plot_glyph(self, **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = 'glyph'
        if 'points' not in kwargs:
            kwargs['points'] = self.mesh
        if 'radius' not in kwargs:
            kwargs['radius'] = self.r * 2.5
        if 'color' not in kwargs:
            kwargs['color'] = self.vertex_color
        if kwargs['name'] not in self.sources or not self.updating:
            P = Points(**kwargs)
            self.add(P)
        else:
            self.update(**kwargs)

    def plot_vectors(self, **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = 'mesh-vectors'
        if 'anchor' not in kwargs:
            kwargs['anchor'] = self.mesh
        if 'color' not in kwargs:
            kwargs['color'] = self.vertex_color
        if kwargs['name'] not in self.sources or not self.updating:
            P = Vectors(**kwargs)
            self.add(P)
            if 'scale_factor' not in kwargs:
                kwargs['scale_factor'] = self.v
                self.update(**kwargs)
        else:
            if 'scale_factor' not in kwargs:
                kwargs['scale_factor'] = self.v
            self.update(**kwargs)

    def hide_edges(self):
        self.hide('edges')

    def remove_edges(self):
        self.remove('edges')

    def hide_faces(self):
        self.hide('faces')

    def remove_faces(self):
        self.remove('faces')
        for name in self._face_sources:
            self.remove(name)
        self._face_sources = []

    def hide_vertices(self):
        self.hide('vertices')

    def remove_vertices(self):
        self.remove('vertices')

    # -------------------------------------------------------------------------
    #                             Virtual plot
    # -------------------------------------------------------------------------

    def selection_on(self):
        M = Edges(self.mesh,
                  name = 'virtual-mesh',
                  color = (0.5,0.5,0.5),
                  opacity = 0.001,)
        self.add(M, pickable=True)
        return self.mesh.E

    def selection_off(self):
        self.remove('virtual-mesh')

    def virtual_edges_on(self):
        if 'virtual-edges' not in self.sources or not self.updating:
            M = Edges(self.mesh,
                      name = 'virtual-edges',
                      color = (0.5,0.5,0.5),
                      opacity = 0.001,)
            self.add(M, pickable=True)
        else:
            self.update('virtual-edges')

    def virtual_edges_off(self):
        self.hide('virtual-edges')

    def virtual_faces_on(self):
        if 'virtual-faces' not in self.sources or not self.updating:
            M = Faces(self.mesh,
                      name = 'virtual-faces',
                      color = (0.5,0.5,0.5),
                      opacity = 0.001,)
            self.add(M, pickable=True)
        else:
            self.update('virtual-faces')

    def virtual_faces_off(self):
        self.hide('virtual-faces')

    def virtual_vertices_on(self):
        if 'virtual-vertices' not in self.sources or not self.updating:
            opacity = 0.001
            if self._show_virtual:
                opacity = 1
            M = Points(self.mesh,
                       name = 'virtual-vertices',
                       glyph_type = 'cube',
                       radius = self.r * 1.1,
                       opacity = opacity)
            self.add(M, pickable=True)
        else:
            self.update('virtual-vertices',
                        radius=self.r*1.1)

    def virtual_vertices_off(self):
        self.hide('virtual-vertices')

    #--------------------------------------------------------------------------
    #                            Selection plot
    #--------------------------------------------------------------------------

    def highlight(self):
        tube_radius = self.g * 1.03
        M = Edges(self.mesh,
                  tube_radius = tube_radius,
                  line_width = 5 * self.scale,
                  color = self.selection_color,
                  shading = False,
                  name = 'highlight')
        self.add(M)

    def highlight_off(self):
        self.remove('highlight')

    def plot_selected_vertices(self):
        if len(self.selected_vertices) > 0:
            if 'selected-vertices' not in self.sources or not self.updating:
                if self.view_mode == 'wireframe' or self.view_mode == '3d':
                    glyph = 'wirefrane'
                    radius = None
                elif self.view_mode == 'solid':
                    radius = self.g * 2.52
                    glyph = 'sphere'
                P = Points(self.mesh,
                           vertex_indices = self.selected_vertices,
                           color = self.selection_color,
                           radius = radius,
                           shading = False,
                           glyph_type = glyph,
                           name = 'selected-vertices',
                           resolution = 18)
                self.add(P)
            else:
                self.update('selected-vertices',
                            vertex_indices=self.selected_vertices,
                            radius = self.g * 2.52)
        else:
            self.hide('selected-vertices')

    def plot_selected_edges(self):
        d = np.zeros(self.mesh.E)
        d[self.selected_edges] = 1
        if self.view_mode == 'wireframe' or self.view_mode == '3d':
            tube_radius = None
        elif self.view_mode == 'solid':
            tube_radius = self.g * 1.03
        if 'selected-edges' not in self.sources or not self.updating:
            M = Edges(self.mesh,
                      edge_data = d,
                      tube_radius = tube_radius,
                      line_width = 5 * self.scale,
                      color = [self.edge_color, self.selection_color],
                      opacity = [0,1],
                      lut_range = '0:+',
                      shading = False,
                      name = 'selected-edges')
            self.add(M)
        else:
            self.update('selected-edges', edge_data=d, tube_radius=tube_radius)

    def plot_selected_faces(self):
        d = np.zeros(self.mesh.F)
        d[self.selected_faces] = 1
        if True:#'selected-faces' not in self.sources or not self.updating:
            M = Faces(self.mesh,
                      face_data = d,
                      color = [self.face_color, self.selection_color],
                      opacity = [0,1],
                      lut_range = '0:+',
                      shading = False,
                      name = 'selected-faces')
            self.add(M)
        else:
            self.update('selected-faces', edge_data=d)

    def hide_selected_vertices(self):
        self.hide('selected-vertices')

    def hide_selected_edges(self):
        self.hide('selected-edges')

    def hide_selected_faces(self):
        self.hide('selected-faces')

    # -------------------------------------------------------------------------
    #                               Selection
    # -------------------------------------------------------------------------

    def select_vertices(self):
        self.select_off()
        self.virtual_vertices_on()
        self.selected_vertices = []

        def callback(p_id):
            if p_id == -1:
                return
            v = p_id // 6
            if v not in self.selected_vertices:
                self.selected_vertices.append(v)
            else:
                self.selected_vertices.remove(v)
            self.plot_selected_vertices()

        self.picker_callback(callback,mode='cell', name='virtual-vertices')

    def on_vertex_selection(self, vertex_callback):
        self.select_off()
        self.virtual_vertices_on()
        self.selected_vertices = []

        def callback(p_id):
            if p_id == -1:
                return
            v = p_id // 6
            self.selected_vertices = [v]
            if vertex_callback is not None:
                vertex_callback(v)
            self.virtual_vertices_on()

        self.picker_callback(callback, mode='cell', name='virtual-vertices')

    def select_vertices_off(self):
        self.virtual_vertices_off()
        self.hide_selected_vertices()
        self.picker_off()
        self.selected_vertices = []

    def select_boundary_vertices(self):
        self.select_off()
        self.virtual_vertices_on()
        self.selected_vertices = []

        def v_callback(p_id):
            if p_id == -1:
                return
            v = p_id // 6
            corners = self.mesh.mesh_corners()
            boundaries = self.mesh.boundary_curves(corner_split=True)
            if v not in corners:
                if v not in self.selected_vertices:
                    for boundary in boundaries:
                        if int(v) in boundary:
                            self.selected_vertices.extend(boundary)
                    self.plot_selected_vertices()
                elif v in self.selected_vertices:
                    for boundary in boundaries:
                        if int(v) in boundary:
                            for w in boundary:
                                self.selected_vertices.remove(w)
                    self.plot_selected_vertices()

        self.picker_callback(v_callback,mode='cell',name='virtual-vertices')

    def select_boundary_vertices_off(self):
        self.select_vertices_off()

    def select_edges(self):
        self.select_off()
        self.virtual_edges_on()
        self.selected_edges = []

        def e_callback(cell_id):
            if cell_id == -1:
                return
            e = cell_id
            if e not in self.selected_edges:
                self.selected_edges.append(e)
            else:
                self.selected_edges.remove(e)
            self.plot_selected_edges()

        self.picker_callback(e_callback, mode='cell', name='virtual-edges')

    def on_edge_selection(self, edge_callback):
        self.select_off()
        self.virtual_edges_on()
        self.selected_edges = []

        def e_callback(cell_id):
            if cell_id == -1:
                return
            e = cell_id
            self.selected_edges = [e]
            if edge_callback is not None:
                edge_callback(e)
            self.virtual_edges_on()

        self.picker_callback(e_callback, mode='cell', name='virtual-edges')

    def select_edges_off(self):
        self.virtual_edges_off()
        self.hide_selected_edges()
        self.picker_off()
        self.selected_edges = []

    def select_faces(self):
        self.select_off()
        self.virtual_faces_on()
        self.selected_faces = []

        def f_callback(cell_id):
            if cell_id == -1:
                return
            f = cell_id
            if f not in self.selected_faces:
                self.selected_faces.append(f)
            else:
                self.selected_faces.remove(f)
            self.plot_selected_faces()

        self.picker_callback(f_callback, mode='cell', name='virtual-faces')

    def on_face_selection(self, face_callback):
        self.virtual_faces_on()
        self.selected_faces = []

        def f_callback(cell_id):
            if cell_id == -1:
                return
            f = cell_id
            self.selected_faces = [f]
            if face_callback is not None:
                face_callback(f)
            self.virtual_faces_on()

        self.picker_callback(f_callback, mode='cell', name='virtual-faces')

    def select_faces_off(self):
        self.virtual_faces_off()
        self.hide_selected_faces()
        self.picker_off()
        self.selected_faces = []

    def clear_selection(self):
        self.hide_selected_vertices()
        self.hide_selected_edges()
        self.hide_selected_faces()
        self.selected_vertices = []
        self.selected_edges = []
        self.selected_faces = []
        self.virtual_edges_off()
        self.virtual_vertices_off()
        self.remove_widgets()

    def select_off(self):
        self.select_edges_off()
        self.select_vertices_off()
        self.select_faces_off()
        self.move_vertex_off()
        self.move_vertices_off()

    #--------------------------------------------------------------------------
    #                               Moving
    #--------------------------------------------------------------------------

    def move_vertex(self, vertex_index, interaction_callback=None,
                   end_callback=None):
        self.remove_widgets()
        S = tvtk.SphereWidget()
        point = self.mesh.vertices[vertex_index]
        S.center = point
        S.radius = self.g * 2
        S.representation = 'surface'
        P = Points(points = point,
                   radius = self.g * 2.52,
                   color = self.selection_color,
                   shading = False,
                   name = 'widget-point')
        self.add(P)
        self.virtual_vertices_on()

        def i_callback(obj, event):
            self.virtual_vertices_off()
            c = obj.GetCenter()
            center = np.array([c[0],c[1],c[2]])
            self.mesh.vertices[vertex_index,:] = center
            self.update('widget-point', points = center)
            if interaction_callback is not None:
                interaction_callback()

        S.add_observer("InteractionEvent", i_callback)

        def e_callback(obj, event):
            if end_callback is not None:
                end_callback()
            self.virtual_vertices_on()
            point = self.mesh.vertices[vertex_index]
            S.center = point
            self.update('widget-point', points = point)

        S.add_observer("EndInteractionEvent", e_callback)
        self.add_widget(S, name='vertex_handler')

    def move_vertex_off(self):
        self.remove_widgets()
        self.remove('widget-point')
        self.virtual_vertices_off()

    def move_vertices(self, interaction_callback=None, start_callback=None,
                      end_callback=None):
        self.select_off()
        self.virtual_vertices_on()
        self.virtual_edges_off()

        def v_callback(p_id):
            self.hide_selected_vertices()
            if p_id == -1:
                return
            v = p_id//6
            self.selected_vertices = [v]
            if start_callback is not None:
                try:
                    start_callback()
                except:
                    pass
            self.move_vertex(v, interaction_callback, end_callback)

        self.picker_callback(v_callback, mode='cell', name='virtual-vertices')

    def move_vertices_off(self):
        self.selected_vertices = []
        self.move_vertex_off()
        self.picker_off()
