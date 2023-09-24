#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

from mayavi.sources.vtk_data_source import VTKDataSource

from mayavi.modules.surface import Surface

from mayavi.filters.tube import Tube

import numpy as np

#------------------------------------------------------------------------------

from geometrylab.vtkplot import plotutilities

# -----------------------------------------------------------------------------

'''polylinesource.py: The plot source for Polylines and Circles'''

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Polyline
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

class Polyline(object):

    def __init__(self, polyline, **kwargs):

        self._polyline      = polyline

        self._edge_data      = kwargs.get('edge_data',      None)

        self._vertex_data    = kwargs.get('vertex_data',    None)

        self._color          = kwargs.get('color',          'white')

        self._opacity        = kwargs.get('opacity',        1)

        self._line_width     = kwargs.get('line_width',     1)

        self._tube_radius    = kwargs.get('tube_radius',    None)

        self._tube_sides     = kwargs.get('tube_sides',     10)

        self._lut_range      = kwargs.get('lut_range',      '-:0:+')

        self._lut_expansion  = kwargs.get('lut_expansion',  2)

        self._reverse_lut    = kwargs.get('reverse_lut',    False)

        self._shading        = kwargs.get('shading',        True)

        self._glossy         = kwargs.get('glossy',         0.3)

        self._sampling       = kwargs.get('sampling',       50)

        self.name            = kwargs.get('name',           'edges')

        self._data_range     = None



    #--------------------------------------------------------------------------
    #                                 Pipeline
    #--------------------------------------------------------------------------
    # - .source  =   VTKDataSource
    # - .data    =     \_Poly_data
    # - .tube    =       None  \_Tube
    # - .module  =        \_______\_Module_manager
    # - .surface =                   \__Surface
    #--------------------------------------------------------------------------

        self._make_data()

        self._make_tube()

        self._module = plotutilities.make_module(self._color,
                                                self._opacity,
                                                self._lut_range,
                                                self._lut_expansion,
                                                self._reverse_lut,
                                                self._data_range)
        self._make_surface()

        self._assemble_pipeline()

        self.on_scene = False

    @property
    def type(self):
        return 'Polyline-source'

    #--------------------------------------------------------------------------
    #                             Data Structure
    #--------------------------------------------------------------------------

    def _make_data(self, **kwargs):
        self._polyline = kwargs.get('polyline', self._polyline) #Hui add
        if self._polyline.type == 'Circle':
            self._polyline.sampling = self._sampling
            self._polyline.make_vertices()
        points = self._polyline.vertices
        try:
            cells = self._polyline.cell_array() #Hui: relate with definition of cell_array() in polyline.py
        except:
            cells = self._polyline.cell_array
        N = cells.shape[0]
        cell_array = tvtk.CellArray()
        cell_array.set_cells(N,cells)
        self._data = tvtk.PolyData(points=points, lines=cell_array)
        if self._edge_data is not None:
            scalars = np.array(self._edge_data)
            self._data.cell_data.scalars = scalars
            self._data.cell_data.scalars.name = self.name
            self._data_range = [np.min(scalars),np.max(scalars)]
        elif self._vertex_data is not None:
            scalars = np.array(self._vertex_data)
            self._data.point_data.scalars = scalars
            self._data.point_data.scalars.name = self.name
            self._data_range = [np.min(scalars),np.max(scalars)]
        else:
            scalars = np.zeros([self._polyline.E])
            self._data.cell_data.scalars = scalars
            self._data.cell_data.scalars.name = self.name
            self._data_range = [0,0]

    #--------------------------------------------------------------------------
    #                                Surface
    #--------------------------------------------------------------------------

    def _make_surface(self):
        self._surface = Surface()
        if self._line_width != None:
            self._surface.actor.property.line_width = self._line_width
        if type(self._opacity) == int or type(self._opacity) == float:
            self._surface.actor.property.opacity = self._opacity
        if not self._shading:
            self._surface.actor.actor.property.lighting = False
        if self._glossy!= 0:
            self._surface.actor.actor.property.specular = 0.7 * self._glossy
            self._surface.actor.actor.property.specular_power = 11*self._glossy

    #--------------------------------------------------------------------------
    #                                 Tube
    #--------------------------------------------------------------------------

    def _make_tube(self):
        self._tube = None
        if self._tube_radius == 'adaptive':
            self._surface.actor.actor.property.representation = 'wireframe'
            self._surface.actor.actor.property.render_lines_as_tubes = True
        elif self._tube_radius is not None:
            self._tube = Tube()
            self._tube.filter.radius = self._tube_radius
            self._tube.filter.number_of_sides = self._tube_sides
            self._tube

    #--------------------------------------------------------------------------
    #                               Pipeline
    #--------------------------------------------------------------------------

    def _assemble_pipeline(self):
        src = VTKDataSource(data=self._data)
        self._module.add_child(self._surface)
        if self._tube is None:
            src.add_child(self._module)
        else:
            self._tube.add_child(self._module)
            src.add_child(self._tube)
        self.source = src

    # -------------------------------------------------------------------------
    #                                 Update
    # -------------------------------------------------------------------------

    def _update_data(self, **kwargs):
        if True:
            self._polyline = kwargs.get('polyline', self._polyline)
            if self._polyline.type == 'Circle':
                self._polyline.sampling = self._sampling
                self._polyline.make_vertices()
            points = self._polyline.vertices
            cells = self._polyline.cell_array()
            N = cells.shape[0]
            cell_array = tvtk.CellArray()
            cell_array.set_cells(N, cells)
            cell_array = tvtk.CellArray()
            cell_array.set_cells(self._polyline.E, cells)
            self._data.set(lines=cell_array)
            self._data.set(points=points)
        self._edge_data = kwargs.get('edge_data', None)
        self._vertex_data = kwargs.get('vertex_data', None)
        if self._edge_data is not None:
            scalars = np.array(self._edge_data)
            self._data.cell_data.scalars = scalars
            self._data.cell_data.scalars.name = self.name
            self._data_range = [np.min(scalars),np.max(scalars)]
        elif self._vertex_data is not None:
            scalars = np.array(self._vertex_data)
            self._data.point_data.scalars = scalars
            self._data.point_data.scalars.name = self.name
            self._data_range = [np.min(scalars),np.max(scalars)]
        else:
            scalars = np.zeros([self._polyline.E])
            self._data.cell_data.scalars = scalars
            self._data.cell_data.scalars.name = self.name
            self._data_range = [0,0]
        self._data.modified()

    def _update_lut_range(self, **kwargs):
        self._lut_range = kwargs.get('lut_range', self._lut_range)
        lut_range = plotutilities.make_lut_range(self._lut_range, self._data_range)
        self._module.scalar_lut_manager.data_range = lut_range

    def _update_color(self, **kwargs):
        if 'color' in kwargs:
            if kwargs['color'] != self._color:
                self._color = kwargs['color']
                plotutilities.make_lut_table(self._module,
                                             self._color,
                                             self._opacity,
                                             self._lut_expansion,
                                             self._reverse_lut)

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

    def update(self, **kwargs):
        self._update_data(**kwargs)
        self._update_lut_range(**kwargs)
        self._update_color(**kwargs)
        self.source.update()

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
