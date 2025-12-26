#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

from mayavi.sources.vtk_data_source import VTKDataSource

from mayavi.modules.glyph import Glyph

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.vtkplot import plotutilities

# -----------------------------------------------------------------------------

'''pointsource.py: Point plot source class, for meshes and arrays'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                   Points
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

'''
-------------------------------------------------------------------------------
Points:

    A point plot object for Mayavi pipeline.

    parameters:

        points:
            a (n,3) 'numpy array' / a 'Mesh'

        point_data:
            a (n) 'numpy array' of points data

        indices:
            a 'numpy array' of indices referred to the Mesh vertex index

        color:
            a (r,g,b) color / a keyword color / a keyword lut-table / a list of
            (r,g,b) and keyword colors. r, g, b in [0,1]

        opacity:
            a [0,1] float / a list of [0,1] floats of the same length of colors

        lut_range:
            '-:0:+'   = symmetric lut range
            '0:+'     = positive lut range
            '-:0'     = negative lut range
            '-:+'     = min-max lut range
            [min,max] = defined lut range

        glyph_type:
            'sphere' / 'cube' / 'axes' / 'cone' / 'wireframe'

-------------------------------------------------------------------------------
'''

class Points(object):

    def __init__(self, points, **kwargs):
        plotutilities.check_arguments(**kwargs)

        self._points         = points

        self._geometry       = None

        self._data_range      = None

        self._vertex_data     = kwargs.get('vertex_data',    None)

        self._indices         = kwargs.get('vertex_indices', None)

        self._color          = kwargs.get('color',          'white')

        self._opacity        = kwargs.get('opacity',        1)

        self._line_width      = kwargs.get('line_width',     1)

        self._scale_factor    = kwargs.get('scale_factor',   1)

        self._glyph_type      = kwargs.get('glyph_type',     'sphere')

        self._radius          = kwargs.get('radius',         0.1)

        self._height          = kwargs.get('height',         0.1)

        self._resolution      = kwargs.get('resolution',     20)

        self._lut_range       = kwargs.get('lut_range',      '-:0:+')

        self._lut_expansion   = kwargs.get('lut_expansion',  2)

        self._reverse_lut     = kwargs.get('reverse_lut',    False)

        self._shading         = kwargs.get('shading',        True)

        self._glossy          = kwargs.get('glossy',         0.3)

        self._scaling         = kwargs.get('scaling',        False)

        self.name             = kwargs.get('name',           'points')

    #--------------------------------------------------------------------------
    #                                 Pipeline
    #--------------------------------------------------------------------------
    # - .source =   VTKDataSource
    # - .data   =     \_Unstructured_grid
    # - .module =        \_Module_manager
    # - .glyph  =           \__Glyph
    #--------------------------------------------------------------------------

        self._make_data()

        self._module = plotutilities.make_module(self._color,
                                                self._opacity,
                                                self._lut_range,
                                                self._lut_expansion,
                                                self._reverse_lut,
                                                self._data_range)

        self._make_glyph()

        self._assemble_pipeline()

        self.on_scene = False

    @property
    def type(self):
        return 'Points-source'

    #--------------------------------------------------------------------------
    #                               Data Structure
    #--------------------------------------------------------------------------

    def _make_data(self):
        if hasattr(self._points, 'type'):
            if self._points.type is 'Mesh' or self._points.type is 'Points':
                self._geometry = self._points
                if self._indices is None:
                    V = self._geometry.V
                    points = self._geometry.vertices
                else:
                    V = len(self._indices)
                    self._indices = np.array(self._indices, dtype=int)
                    points = self._geometry.vertices[self._indices,:]
        else:
            try:
                self._points.shape
            except AttributeError:
                self._points = np.array(self._points)
            if len(self._points.shape) == 1:
                points = np.array([self._points])
            else:
                points = self._points
        self._data = tvtk.PolyData(points=points)
        if self._vertex_data is None:
            V = points.shape[0]
            self._data.point_data.scalars = np.zeros(V)
            self._data_range = [0,0]
        else:
            scalars = np.array(self._vertex_data)
            self._data.point_data.scalars = scalars
            self._data_range = [np.min(scalars),np.max(scalars)]
        self._data.point_data.scalars.name = self.name

    #--------------------------------------------------------------------------
    #                                Glyph
    #--------------------------------------------------------------------------

    def _make_glyph(self):
        self._glyph = Glyph()
        if self._scaling:
            self._glyph.glyph.scale_mode = 'scale_by_scalar'
        else:
            self._glyph.glyph.scale_mode = 'data_scaling_off'
        self._glyph.actor.property.opacity = self._opacity
        #print(glyph.glyph.glyph_source.glyph_dict)

        if self._glyph_type == 'wireframe':
            self._glyph.actor.actor.property.representation = 'points'
            self._glyph.actor.actor.property.point_size = 1
            self._glyph.actor.actor.property.render_points_as_spheres = True
            g_src = self._glyph.glyph.glyph_source.glyph_dict['glyph_source2d']
            self._glyph.glyph.glyph_source.glyph_source.dash = True
            self._glyph.glyph.glyph_source.glyph_source.glyph_type = 'none'
            self._glyph.glyph.glyph_source.glyph_source.scale = 0
            self._glyph.actor.property.line_width = self._line_width

        elif self._glyph_type == 'cube':
            self._glyph.actor.actor.property.representation = 'surface'
            g_src = self._glyph.glyph.glyph_source.glyph_dict['cube_source']
            g_src.x_length = 2*self._radius
            g_src.y_length = 2*self._radius
            g_src.z_length = 2*self._radius

        elif self._glyph_type == 'sphere':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['sphere_source']
            g_src.radius = self._radius
            g_src.phi_resolution = self._resolution
            g_src.theta_resolution = self._resolution

        elif self._glyph_type == 'cone':
            self._glyph.actor.actor.property.representation = 'surface'
            g_src = self._glyph.glyph.glyph_source.glyph_dict['cone_source']
            g_src.resolution = self._resolution
            g_src.radius = self._radius
            g_src.height = self._height
            g_src.direction =  np.array([0, 0, 1])
            g_src.center = np.array([0, 0, -self._height/2])


        elif self._glyph_type == 'axes':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['axes']
            self._glyph.actor.property.line_width = self._line_width

        self._glyph.glyph.glyph_source.glyph_source = g_src
        self._glyph.glyph.glyph.scale_factor = self._scale_factor

        if not self._shading:
            self._glyph.actor.actor.property.lighting = False

        self._glyph.actor.actor.property.specular = min(self._glossy, 1)
        sp = 21. * self._glossy + 0.001
        self._glyph.actor.actor.property.specular_power = sp

    #--------------------------------------------------------------------------
    #                                Pipeline
    #--------------------------------------------------------------------------

    def _assemble_pipeline(self):
        src = VTKDataSource(data=self._data)
        self._module.add_child(self._glyph)
        src.add_child(self._module)
        self.source = src

    # -------------------------------------------------------------------------
    #                                 Update
    # -------------------------------------------------------------------------

    def _update_data(self, **kwargs):
        self._indices = kwargs.get('vertex_indices', self._indices)
        if self._geometry is not None:
            if self._points.type is 'Mesh' or self._points.type is 'Points':
                if self._indices is not None:
                    if type(self._indices) == list:
                        self._indices = np.array(self._indices)
                    points = self._geometry.vertices[self._indices,:]
                elif 'points' in kwargs:
                    points = kwargs['points']
                    if type(points) == list:
                        points = np.array(points)
                    if len(points.shape) == 1:
                        points = np.array([points])
                else:
                    points = self._geometry.vertices
        elif 'points' in kwargs:
            points = kwargs['points']
            if type(points) == list:
                points = np.array(points)
            if len(points.shape) == 1:
                points = np.array([points])
        else:
            points = self._points
        self._data.points=points ##Hui replace .set to .points
        self._vertex_data = kwargs.get('vertex_data', self._vertex_data)
        if self._vertex_data is not None:
            scalars = np.array(self._vertex_data)
            self._data.point_data.scalars = scalars
            self._data.point_data.scalars.name = self.name
            self._data_range = [np.min(scalars), np.max(scalars)]
        else:
            V = points.shape[0]
            self._data.point_data.scalars = np.zeros(V)
            self._data.point_data.scalars.name = self.name
            self._data_range = [0,0]
        self._data.modified()

    def _update_lut_range(self, **kwargs):
        self._lut_range = kwargs.get('lut_range', self._lut_range)
        lut_range = plotutilities.make_lut_range(self._lut_range, self._data_range)
        self._module.scalar_lut_manager.data_range = lut_range

    def _update_glossy(self, **kwargs):
        self._glossy = kwargs.get('glossy', self._glossy)
        self._glyph.actor.actor.property.specular = min(self._glossy, 1)
        sp = 21. * self._glossy + 0.001
        self._glyph.actor.actor.property.specular_power = sp

    def _update_radius(self, **kwargs):
        g_src = self._glyph.glyph.glyph_source.glyph_source
        if self._glyph_type is 'cube':
            g_src.x_length = 2*self._radius
            g_src.y_length = 2*self._radius
            g_src.z_length = 2*self._radius
        elif self._glyph_type is 'sphere':
            self._radius = kwargs.get('radius', self._radius)
            g_src.radius = self._radius

    def _update_color(self, **kwargs):
        if 'color' in kwargs:
            if kwargs['color'] != self._color:
                self._color = kwargs['color']
                plotutilities.make_lut_table(self._module,
                                           self._color,
                                           self._opacity,
                                           self._lut_expansion,
                                           self._reverse_lut)


    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    def update(self, **kwargs):
        self._update_data(**kwargs)
        self._update_lut_range(**kwargs)
        self._update_color(**kwargs)
        self._update_glossy(**kwargs)
        self._update_radius(**kwargs)
        # --------------------------- bug fix ---------------------------------
        # plot 2 times for updation of cone center
        if self._glyph_type == 'cone':
            g_src = self._glyph.glyph.glyph_source.glyph_source
            g_src.center = np.array([0, 0, -self._height/2])
        # ---------------------------------------------------------------------
        self.source.update()

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    #                                Legend  ## Hui add
    # -------------------------------------------------------------------------

    def show_legend(self, **kwargs):
        label = kwargs.get('label', 'vertices')
        self._module.scalar_lut_manager.show_legend = True
        self._module.scalar_lut_manager.data_name = label
        #print(dir(self._module.scalar_lut_manager.scalar_bar))

    def hide_legend(self):
        self._module.scalar_lut_manager.show_legend = False

