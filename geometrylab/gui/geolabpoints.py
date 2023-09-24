#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

import numpy as np

from pyface.image_resource import ImageResource

from traits.api import HasTraits, Instance, Property, Enum, Button,String,\
                       on_trait_change, Float, Bool, Int,Constant, ReadOnly,\
                       List, Array, Range, Str, Color

from traitsui.api import View, Item, HSplit, VSplit, InstanceEditor, HGroup,\
                         Group, ListEditor, Tabbed, VGroup, CheckListEditor,\
                         ArrayEditor, Action, ToolBar, Separator,EnumEditor,\
                         ListStrEditor, ColorEditor, Controller

#------------------------------------------------------------------------------

from geometrylab.vtkplot.pointsplotmanager import PointsPlotManager

#------------------------------------------------------------------------------

'''-'''

__author__ = 'Davide Pellis'

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Mesh PlotManager
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class GPHandler(Controller):

    def close(self, info, is_ok):
        info.object._closed = True
        Controller.close(self ,info, is_ok)
        return True


class GeolabPoints(PointsPlotManager):

    _handler = GPHandler()

    _closed = Bool(True)

    x = Int(1)

    scale_range = Range(-10, 10, 0, label = 'scale')

    vector_scale_range = Range(-10, 10, 0, label = 'vectors_scale')

    vertices_glossy_range = Range(0., 1., 0.5)

    vertices_color_select = Color((122, 163, 230), label='vertices color')

    vertices_callbacks = ['color', 'none']

    vertices_plot = Str('color', label = 'Vertices')

    vertices_data_range = Range(0.00, 1.00, 0.20)

    view = View(
                VGroup(VGroup(
                       Item('scale_range',
                            resizable=True),
                       Item('vector_scale_range',
                            resizable=True),
                       show_border=True),
                       HGroup(VGroup(Item('vertices_plot',
                                          resizable=True,
                                          editor=CheckListEditor(values=
                                                           vertices_callbacks),),
                                     ),
                              VGroup(Item('vertices_color_select',
                                          show_label=False,
                                          resizable=False,
                                          enabled_when='vertices_plot == "color"',
                                          editor=ColorEditor(),),
                                     ),
                               ),
                       HGroup(Item('vertices_glossy_range',
                                   resizable=True),
                               ),
                       HGroup(Item('vertices_data_range',
                                   resizable=True),
                              enabled_when=('face_plot == "gaussian curv." \
                                             or face_plot == "mean curv."'),
                               ),
                       Item('_'),
                       show_border=True),
                handler = _handler,
                title= 'Points plot settings',
                icon = ImageResource('img/new/plotmanager.png'),)

    def __init__(self, *args):
        PointsPlotManager.__init__(self, *args)

        self._vertices_callbacks = {'color': self.plot_plain_vertices,
                                    'none': self.hide_vertices}

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

    #--------------------------------------------------------------------------
    #                                    Plot
    #--------------------------------------------------------------------------

    def add_vertices_callback(self, callback, name):
        self._vertices_callbacks[name] = callback
        if name not in self.vertices_callbacks:
            self.vertices_callbacks.append(name)

    @on_trait_change('point_plot')
    def point_plot_change(self):
        self.remove_vertices()
        self.update_plot()

    def update_plot(self):
        PointsPlotManager.update_plot(self)
        vertices_plot = self._vertices_callbacks[self.vertices_plot]
        vertices_plot()

    #--------------------------------------------------------------------------
    #                             Predefined Plots
    #--------------------------------------------------------------------------

    def plot_plain_vertices(self):
        self.plot_vertices(glossy=self.vertices_glossy_range)

    #--------------------------------------------------------------------------
    #                                 Select
    #--------------------------------------------------------------------------

    def select_off(self):
        PointsPlotManager.select_off(self)

    #--------------------------------------------------------------------------
    #                                 Scales
    #--------------------------------------------------------------------------

    @on_trait_change('vertices_data_range, smooth_vertex_data, \
                     vertices_glossy_range')
    def plot_updated(self):
        self.update_plot()

    @on_trait_change('vertices_plot')
    def vertices_plot_change(self):
        self.remove_vertices()
        self.clear() #uneff
        self.update_plot()

    @on_trait_change('scale_range')
    def set_scale(self):
        self.scale = 4**(0.1 * self.scale_range)
        self.glyph_scale = 4**(0.1 * self.scale_range)
        self.update_plot()

    @on_trait_change('vector_scale_range')
    def set_vector_scale(self):
        scale = 4**(0.2 * self.vector_scale_range)
        self.vector_scale = scale
        self.update_plot()

    @on_trait_change('vertices_color_select')
    def set_color(self):
        vc = self.vertices_color_select.getRgb()
        self.vertex_color = (vc[0], vc[1], vc[2])
        self.update_plot()
