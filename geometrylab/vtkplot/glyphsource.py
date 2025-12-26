#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

from mayavi.filters.poly_data_normals import PolyDataNormals

from mayavi.sources.vtk_data_source import VTKDataSource

from mayavi.modules.surface import Surface

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.vtkplot import plotutilities

from geometrylab.vtkplot import glyphs

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                  Length Vector
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Glyph3D:

    def __init__(self, points, **kwargs):

        self.points         = points

        self.geometry       = None

        self.data_range      = None

        self.vectors        = kwargs.get('vectors',       None)

        self.scalars        = kwargs.get('scalars',       None)

        self.glyph_type     = kwargs.get('glyph_type',    '3D-arrow')

        self.anchor_mode    = kwargs.get('anchor_mode',   'vertex')

        self.offset         = kwargs.get('offset',        0)

        self.color          = kwargs.get('color',         'blue-red')

        self.opacity        = kwargs.get('opacity',       1)

        self.glyph_type     = kwargs.get('glyph_type',    '3D-arrow')

        self.position       = kwargs.get('position',      'head')

        self.scale_factor   = kwargs.get('scale_factor',   1)

        self.radius         = kwargs.get('radius',         0.1)

        self.tip_radius     = kwargs.get('tip_radius',     None)

        self.tip_length     = kwargs.get('tip_length',     None)

        self.shaft_radius   = kwargs.get('shaft_radius',   None)

        self.resolution     = kwargs.get('resolution',     20)

        self.lut_range      = kwargs.get('lut_range',      '-:+')

        self.lut_expansion  = kwargs.get('lut_expansion',  2)

        self.reverse_lut    = kwargs.get('reverse_lut',    False)

        self.shading        = kwargs.get('shading',        True)

        self.smooth         = kwargs.get('smooth',         True)

        self.glossy         = kwargs.get('glossy',         0.3)

        self.name           = kwargs.get('name',           'vectors')

    #--------------------------------------------------------------------------
    #                                 Pipeline
    #--------------------------------------------------------------------------
    # - .source  =   VTKDataSource
    # - .data    =     \_Unstructured_grid
    # - .normals =       None  \_Poly_data_normals
    # - .module  =        \_______\_Module_manager
    # - .surface =                   \__Surface
    #--------------------------------------------------------------------------

        self._make_data()

        self._make_normals()

        self.module = plotutilities.make_module(self.color,
                                               self.opacity,
                                               self.lut_range,
                                               self.lut_expansion,
                                               self.reverse_lut,
                                               self.data_range)

        self._make_surface()

        self._assemble_pipeline()

    @property
    def type(self):
        return 'Glyph3D-source'

    #--------------------------------------------------------------------------
    #                               Data Structure
    #--------------------------------------------------------------------------

    def _make_data(self):
        try:
            self.points.shape
        except AttributeError:
            self.points = np.array(self.points)
        if len(self.points.shape) == 1:
            self.points = np.array([self.points])
        if self.glyph_type == '3D-arrow':
            if hasattr(self.points, 'type'):
                if self.points.type == 'mesh':
                    self.geometry = self.points
                    if self.anchor_mode == 'vertex':
                        points = self.geometry.vertices
                    elif self.anchor_mode == 'face':
                        points = self.geometry.face_barycenters()
                    elif self.anchor_mode == 'edge':
                        points = self.geometry.edge_mid_points()
            else:
                points = self.points
            if self.vectors is None:
                self.vectors = np.zeros(self.points.shape)
                self.vectors[:,2] = -1
            mesh = glyphs.Arrows(self.vectors, points,
                                 shaft_radius=self.shaft_radius,
                                 tip_length=self.tip_length,
                                 tip_radius=self.tip_radius,
                                 position=self.position,
                                 offset=self.offset,
                                 resolution=self.resolution,
                                 scale_factor=self.scale_factor)
        elif self.glyph_type == 'sphere':
            mesh = glyphs.Spheres(self.points, self.radius,
                                  resolution=self.resolution)
        else:
            mesh = None
        cell_types = np.array([tvtk.Triangle().cell_type,
                               tvtk.Quad().cell_type,
                               tvtk.Polygon().cell_type])
        cells, types = mesh.cell_arrays()
        cell_array = tvtk.CellArray()
        cell_array.set_cells(mesh.F, cells)
        cell_types = cell_types[types]
        self.data = tvtk.UnstructuredGrid(points = mesh.vertices)
        self.data.set_cells(cell_types, cell_array)
        if self.scalars is not None:
            scalars = np.repeat(self.scalars, mesh.G)
            self.data_range = [np.min(self.scalars), np.max(self.scalars)]
        else:
            if self.glyph_type == '3D-arrow':
                data = mesh.lengths
            elif self.glyph_type == 'sphere':
                data = mesh.radius
            else:
                data = np.zeros(mesh.V)
            scalars = np.repeat(data, mesh.G)
            self.data_range = [np.min(data), np.max(data)]
        self.data.point_data.scalars = scalars
        self.data.point_data.scalars.name = self.name

    #--------------------------------------------------------------------------
    #                                Surface
    #--------------------------------------------------------------------------

    def _make_surface(self):
        self.surface = Surface()
        if not self.shading:
            self.surface.actor.actor.property.lighting = False
        if type(self.opacity) == int or type(self.opacity) == float:
            self.surface.actor.property.opacity = self.opacity
            if self.opacity < 1:
                self.surface.actor.actor.property.specular = 0.8
        if self.glossy != 0:
            self.surface.actor.actor.property.specular = 0.7 * self.glossy
            self.surface.actor.actor.property.specular_power = 11*self.glossy

    #--------------------------------------------------------------------------
    #                                Normals
    #--------------------------------------------------------------------------

    def _make_normals(self):
        if self.smooth:
            self.normals = PolyDataNormals()
        else:
            self.normals = None

    #--------------------------------------------------------------------------
    #                                Pipeline
    #--------------------------------------------------------------------------

    def _assemble_pipeline(self):
        src = VTKDataSource(data=self.data)
        self.module.add_child(self.surface)
        if self.normals is None:
            src.add_child(self.module)
        else:
            self.normals.add_child(self.module)
            src.add_child(self.normals)
        self.source = src
