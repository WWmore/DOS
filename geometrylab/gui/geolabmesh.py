#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from pyface.image_resource import ImageResource

from traits.api import on_trait_change, Bool, Range, Str, Color

from traitsui.api import View, Item, VGroup, HGroup, CheckListEditor, \
    ColorEditor, Controller

# ------------------------------------------------------------------------------

from geometrylab.vtkplot import MeshPlotManager
# ------------------------------------------------------------------------------

'''-'''

__author__ = 'Davide Pellis, Hui Wang'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Mesh PlotManager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class GMHandler(Controller):

    def close(self, info, is_ok):
        info.object._closed = True
        Controller.close(self, info, is_ok)
        return True


class GeolabMesh(MeshPlotManager):
    _handler = GMHandler()

    _closed = Bool(True)

    plot_scale = Range(-10, 10, 0)

    set_vector_scale = Range(-10, 10, 0, label='vector scale')

    face_glossy_range = Range(0., 1., 0.5)

    vertex_color_select = Color((122, 163, 230), width=20, label='vertex_color')

    edge_color_select = Color((122, 163, 230), width=20, label='edge_color')

    face_color_select = Color((122, 163, 230), width=20)

    vertex_plot = Str('none', label='Vertices')

    vertex_callbacks = ['color', 'none']

    edge_plot = Str('color', label='Edges')

    edge_callbacks = ['color', 'wireframe', 'axial force', 'none']

    face_plot = Str('face planarity', label='Faces')

    face_callbacks = ['color', 'face planarity', 'gaussian curv.',
                      'mean curv.', 'none']

    show_supported_vertices = Bool(False, label='Supports')

    show_fixed_vertices = Bool(False, label='Fixed')

    show_loads = Bool(False, label='Loads')

    show_curvature_directions = Bool(False, label='Curv dir')

    face_data_range = Range(0.00, 1.00, 0.20)

    smooth_face_data = Bool(False, label='smooth')

    smooth_faces = Bool(False, label='smooth')

    show_vertex_legend = Bool(False, label='legend')

    show_edge_legend = Bool(False, label='legend')

    show_face_legend = Bool(False, label='legend')

    view = View(
        VGroup(VGroup(
            Item('plot_scale',
                 resizable=True),
            Item('set_vector_scale',
                 resizable=True),
            show_border=True),
            VGroup(
                HGroup(VGroup(Item('vertex_plot',
                                   resizable=True,
                                   editor=CheckListEditor(values=
                                                          vertex_callbacks), ),
                              Item('edge_plot',
                                   resizable=True,
                                   editor=CheckListEditor(values=
                                                          edge_callbacks), ),
                              Item('face_plot',
                                   resizable=True,
                                   editor=CheckListEditor(values=
                                                          face_callbacks), ),
                              ),
                       VGroup(Item('vertex_color_select',
                                   show_label=False,
                                   resizable=False,
                                   enabled_when='vertex_plot == "color"',
                                   editor=ColorEditor()),
                              Item('edge_color_select',
                                   show_label=False,
                                   resizable=False,
                                   enabled_when='edge_plot == "color"',
                                   editor=ColorEditor(), ),
                              Item('face_color_select',
                                   show_label=False,
                                   resizable=False,
                                   enabled_when='face_plot == "color"',
                                   editor=ColorEditor()),
                              ),
                       VGroup('show_vertex_legend',
                              '_',
                              'show_edge_legend',
                              '_',
                              'show_face_legend',
                              ),
                       ),
                show_border=True),
            VGroup(
                HGroup(Item('face_glossy_range',
                            resizable=True),
                       Item('smooth_faces'),
                       enabled_when=('face_plot == "color" \
                                             or face_plot == "face planarity" \
                                             or face_plot == "gaussian curv." \
                                             or face_plot == "mean curv."'),
                       ),
                HGroup(Item('face_data_range',
                            resizable=True),
                       Item('smooth_face_data'),
                       enabled_when=('face_plot == "gaussian curv." \
                                             or face_plot == "mean curv."'),
                       ),
                show_border=True),
            VGroup(
                HGroup(VGroup('show_supported_vertices',
                              ),
                       VGroup('show_fixed_vertices',
                              ),
                       VGroup('show_loads',
                              ),
                       VGroup('show_curvature_directions',
                              ),
                       ),
                show_border=True),
            show_border=True),
        handler=_handler,
        title='Mesh plot settings',
        icon=ImageResource('img/new2/plotmanager.png'), )

    def __init__(self, *args):
        MeshPlotManager.__init__(self, *args)

        self._vertex_callbacks = {'color': self.plot_plain_vertices,
                                  'none': self.hide_vertices}

        self._edge_callbacks = {'color': self.plot_plain_edges,
                                'wireframe': self.plot_wireframe_edges,
                                'axial force': self.plot_axial_force,
                                'none': self.hide_edges}

        self._face_callbacks = {'color': self.plot_plain_faces,
                                'face planarity': self.plot_face_planarity,
                                'gaussian curv.': self.plot_gaussian_curvature,
                                'mean curv.': self.plot_mean_curvature,
                                'none': self.hide_faces}

    def start(self):
        if not self._closed:
            pass
            self._handler.info.ui.dispose()
        self.configure_traits()
        self._closed = False

    def close(self):
        try:
            self._closed = True
            self._handler.info.ui.dispose()
        except:
            pass

    # --------------------------------------------------------------------------
    #                                    Plot
    # --------------------------------------------------------------------------

    def add_vertex_callback(self, callback, name):
        self._vertex_callbacks[name] = callback
        if name not in self.vertex_callbacks:
            self.vertex_callbacks.append(name)

    def add_edge_callback(self, callback, name):
        self._edge_callbacks[name] = callback
        if name not in self.edge_callbacks:
            self.edge_callbacks.append(name)

    def add_face_callback(self, callback, name):
        self._face_callbacks[name] = callback
        if name not in self.face_callbacks:
            self.face_callbacks.append(name)

    @on_trait_change('vertex_plot')
    def vertex_plot_change(self):
        self.remove_vertices()
        self.clear()  # uneff
        self.update_plot()

    @on_trait_change('edge_plot')
    def edge_plot_change(self):
        self.remove_edges()
        self.clear()  # uneff
        self.update_plot()

    @on_trait_change('face_plot, smooth_faces')
    def face_plot_changed(self):
        self.remove_faces()
        self.clear()
        self.update_plot()

    def update_plot(self, **kwargs):
        # self.disable_render()
        MeshPlotManager.update_plot(self, **kwargs)
        self.plot_supported_vertices()
        self.plot_fixed_vertices()
        self.plot_loads()
        self.plot_curvature_directions()
        vertex_plot = self._vertex_callbacks[self.vertex_plot]
        vertex_plot()
        edge_plot = self._edge_callbacks[self.edge_plot]
        edge_plot()
        face_plot = self._face_callbacks[self.face_plot]
        face_plot()
        self.make_face_legend()
        # self.enable_render()

    # --------------------------------------------------------------------------
    #                             Predefined Plots
    # --------------------------------------------------------------------------

    def plot_plain_vertices(self):
        self.plot_vertices()

    def plot_plain_edges(self):
        self.plot_edges()

    def plot_axial_force(self):
        edge_data = -self.mesh.axial_forces
        self.plot_edges(edge_data=edge_data,
                        color='bwr_e',
                        lut_range='-:0:+')

    def plot_wireframe_edges(self):
        self.plot_edges(color='black',
                        tube_radius=None)

    def plot_plain_faces(self):
        self.plot_faces(smooth=self.smooth_faces,
                        glossy=self.face_glossy_range)

    def plot_face_planarity(self):
        P = self.mesh.face_planarity()
        val = np.max(P) * self.face_data_range
        val = 0.02
        self.plot_faces(face_data=P,
                        smooth=self.smooth_faces,
                        glossy=self.face_glossy_range,
                        lut_range=[0, val],
                        color='blue-red')

    def plot_gaussian_curvature(self):
        K = self.mesh.principal_curvatures(area_normalization=True)
        vertex_data = K[0] * K[1]
        if self.smooth_face_data:
            vertex_data = self.mesh.smooth_vertex_data(vertex_data)
        val = np.max(np.abs(vertex_data)) * self.face_data_range
        self.plot_faces(vertex_data=vertex_data,
                        smooth=self.smooth_faces,
                        glossy=self.face_glossy_range,
                        color='blue-red',
                        lut_range=[-val, val])

    def plot_mean_curvature(self):
        K = self.mesh.principal_curvatures(area_normalization=True)
        vertex_data = (K[0] + K[1]) / 2
        if self.smooth_face_data:
            vertex_data = self.mesh.smooth_vertex_data(vertex_data)
        val = np.max(np.abs(vertex_data)) * self.face_data_range
        self.plot_faces(vertex_data=vertex_data,
                        smooth=self.smooth_faces,
                        glossy=self.face_glossy_range,
                        color='blue-red',
                        lut_range=[-val, val])

    @on_trait_change('show_fixed_vertices')
    def plot_fixed_vertices(self):
        v = self.mesh.fixed_vertices
        if self.show_fixed_vertices and len(v) != 0:
            r = self.g
            self.plot_glyph(vertex_indices=v,
                            glyph_type='axes',
                            color='orange',
                            shading=False,
                            line_width=5,
                            scale_factor=r ** 0.5 * 2,
                            name='fixed')
        else:
            self.hide('fixed')

    @on_trait_change('show_supported_vertices')
    def plot_supported_vertices(self):
        v = self.mesh.constrained_vertices
        if self.show_supported_vertices and len(v) != 0:
            r = self.r
            self.plot_glyph(vertex_indices=v,
                            glyph_type='sphere',
                            resolution=20,
                            radius=r * 2,
                            color='m',
                            name='supports')
        else:
            self.hide('supports')

    @on_trait_change('show_loads')
    def plot_loads(self):
        if self.show_loads:
            self.plot_vectors(vectors=self.mesh.loads(),
                              color='blue-red',
                              glyph_type='3D-arrow',
                              position='head',
                              name='loads')
        else:
            self.hide('loads')

    @on_trait_change('show_curvature_directions')
    def plot_curvature_directions(self):
        if self.show_curvature_directions:
            K = self.mesh.principal_curvatures()
            self.plot_vectors(vectors=K[2],
                              position='center',
                              glyph_type='line',
                              line_width=2,
                              color='b',
                              name='c-d1')
            self.plot_vectors(vectors=K[3],
                              position='center',
                              glyph_type='line',
                              line_width=2,
                              color='r',
                              name='c-d2')
        else:
            self.hide('c-d1')
            self.hide('c-d2')

    @on_trait_change('show_face_legend')
    def make_face_legend(self):
        try:
            if self.show_face_legend:
                self.get_source('faces').show_legend(label=self.face_plot)
            else:
                self.get_source('faces').hide_legend()
        except:
            pass

    @on_trait_change('show_vertex_legend') # Hui add
    def make_vertex_legend(self):
        try:
            if self.show_vertex_legend:
                self.get_source('vertices').show_legend(label=self.vertex_plot)
            else:
                self.get_source('vertices').hide_legend()
        except:
            pass
        
    @on_trait_change('show_edge_legend') # Hui add
    def make_edge_legend(self):
        try:
            if self.show_edge_legend:
                self.get_source('edges').show_legend(label=self.edge_plot)
            else:
                self.get_source('edges').hide_legend()
        except:
            pass
    # --------------------------------------------------------------------------
    #                                 Select
    # --------------------------------------------------------------------------

    def select_off(self):
        super(GeolabMesh, self).select_off()

    # --------------------------------------------------------------------------
    #                                 Scales
    # --------------------------------------------------------------------------

    @on_trait_change('face_data_range, smooth_face_data, face_glossy_range')
    def plot_updated(self):
        self.update_plot()

    @on_trait_change('plot_scale')
    def set_scale(self):
        # self.handler.set_state(None)
        self.scale = 4 ** (0.1 * self.plot_scale)
        self.glyph_scale = 4 ** (0.1 * self.plot_scale)
        self.update_plot()

    @on_trait_change('set_vector_scale')
    def set_set_vector_scale(self):
        # self.handler.set_state(None)
        scale = 4 ** (0.2 * self.set_vector_scale)
        self.vector_scale = scale
        self.update_plot()

    @on_trait_change('face_color_select, edge_color_select, vertex_color_select')
    def set_color(self):
        vc = self.vertex_color_select.getRgb()
        ec = self.edge_color_select.getRgb()
        fc = self.face_color_select.getRgb()
        self.vertex_color = (vc[0], vc[1], vc[2])
        self.edge_color = (ec[0], ec[1], ec[2])
        self.face_color = (fc[0], fc[1], fc[2])
        self.update_plot()
