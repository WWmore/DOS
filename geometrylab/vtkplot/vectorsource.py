#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from tvtk.api import tvtk

from mayavi.sources.vtk_data_source import VTKDataSource

from mayavi.modules.vectors import Vectors as vtk_vectors

from mayavi.modules.glyph import Glyph

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.vtkplot import plotutilities

# -----------------------------------------------------------------------------

'''vectorsource.py: The vectors plot source class'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                   Vectors
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class Vectors:

    def __init__(self, vectors, **kwargs):
        plotutilities.check_arguments(**kwargs)

        self._vectors        = vectors

        self._geometry       = None

        self._anchor         = kwargs.get('anchor',    np.zeros(vectors.shape))

        self._scalars        = kwargs.get('scalars',        None)

        self._anchor_mode     = kwargs.get('anchor_mode',    'vertex')

        self._color          = kwargs.get('color',          'blue-red')

        self._opacity        = kwargs.get('opacity',        1)

        self._line_width     = kwargs.get('line_width',     1)

        self._glyph_type     = kwargs.get('glyph_type',     '3D-arrow')

        self._position       = kwargs.get('position',       'head')

        self._scaling        = kwargs.get('scaling',        True)

        self._scale_factor   = kwargs.get('scale_factor',    1)

        self._radius         = kwargs.get('radius',         None)

        self._tip_radius     = kwargs.get('tip_radius',     None)

        self._tip_length     = kwargs.get('tip_length',     None)

        self._shaft_radius   = kwargs.get('shaft_radius',   None)

        self._height         = kwargs.get('height',         None)

        self._resolution     = kwargs.get('resolution',     20)

        self._lut_range      = kwargs.get('lut_range',      '-:+')

        self._lut_expansion  = kwargs.get('lut_expansion',  2)

        self._reverse_lut     = kwargs.get('reverse_lut',    False)

        self._shading        = kwargs.get('shading',        True)

        self._glossy         = kwargs.get('glossy',         0.3)

        self.name            = kwargs.get('name',           'vectors')

    #--------------------------------------------------------------------------
    #                                 Pipeline
    #--------------------------------------------------------------------------
    # - .source =   VTKDataSource
    # - .data   =     \_Unstructured_grid
    # - .module =        \_Module_manager
    # - .glyph  =           \__Glyph
    #--------------------------------------------------------------------------

        self._make_data()

        self._make_glyph()

        self._module = plotutilities.make_vector_module(self._color,
                                                      self._opacity,
                                                      self._lut_range,
                                                      self._lut_expansion,
                                                      self._reverse_lut)

        self._assemble_pipeline()

        self.on_scene = False

        @property
        def type(self):
            return 'Vectors-source'
    #--------------------------------------------------------------------------
    #                               Data Structure
    #--------------------------------------------------------------------------

    def _make_data(self):
        try:
            self._vectors.shape
        except:
            self._vectors = np.array(self._vectors)
        if len(self._vectors.shape) == 1:
            self._vectors = np.array([self._vectors])
        norm = np.linalg.norm(self._vectors, axis=1)
        if hasattr(self._anchor, 'vertices'):
            if self._anchor.type is 'Mesh':
                self._geometry = self._anchor
                if self._anchor_mode == 'vertex':
                    points = self._geometry.vertices
                elif self._anchor_mode == 'face':
                    points = self._geometry.face_barycenters()
                elif self._anchor_mode == 'edge':
                    points = self._geometry.edge_mid_points()
            if self._anchor.type is 'Points':
                self._geometry = self._anchor
                points = self._geometry.vertices
        else:
            try:
                self._anchor.shape
            except:
                self._anchor = np.array(self._anchor)
            if len(self._anchor.shape) == 1:
                self._anchor = np.array([self._anchor])
            points = self._anchor
        mask = np.where(norm > 0)[0]
        points = points[mask]
        vectors = self._vectors[mask]
        self._data = tvtk.UnstructuredGrid(points=points)
        if self._scalars is not None:
            scalars = self._scalars[mask]
            self._data.point_data.scalars = scalars
            self._data.point_data.scalars.name = 'scalars'
        self._data.point_data.vectors = vectors
        self._data.point_data.vectors.name = self.name

    #--------------------------------------------------------------------------
    #                                Glyph
    #--------------------------------------------------------------------------

    def _make_glyph(self):

        self._glyph = Glyph()
        self._glyph.glyph.scale_mode = 'scale_by_vector'
        self._glyph.glyph.color_mode = 'color_by_vector'
        self._glyph.actor.property.opacity = self._opacity
        #print(vector.glyph.glyph_source.glyph_dict)

        if self._glyph_type == '3D-arrow':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['arrow_source']
            g_src.shaft_resolution = self._resolution
            g_src.tip_resolution = self._resolution
            if self._shaft_radius != None:
                g_src.shaft_radius = self._shaft_radius
            if self._tip_length != None:
                g_src.tip_length = self._tip_length
            if self._tip_radius != None:
                g_src.tip_radius = self._tip_radius
            self._glyph.glyph.glyph_source.glyph_source = g_src

        elif self._glyph_type == 'cylinder':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['cylinder_source']
            g_src.resolution = self._resolution
            g_src.radius = self._radius
            g_src.height = self._height
            self._glyph.glyph.glyph_source.glyph_source = g_src

        elif self._glyph_type == 'line':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['glyph_source2d']
            g_src.dash = True
            g_src.glyph_type = 'none'
            g_src.scale = self._scale_factor
            #g_src.line_width = self._line_width
            self._glyph.glyph.glyph_source.glyph_source = g_src

        elif self._glyph_type == 'cone':
            self._glyph.actor.actor.property.representation = 'surface'
            g_src = self._glyph.glyph.glyph_source.glyph_dict['cone_source']
            g_src.resolution = self._resolution
            g_src.radius = self._radius
            g_src.height =self._height
            self._glyph.glyph.glyph_source.glyph_source = g_src

        elif self._glyph_type == 'axes':
            g_src = self._glyph.glyph.glyph_source.glyph_dict['axes']
            g_src.line_width = self._line_width
            #g_src.scale_factor = self._scale_factor
            self._glyph.glyph.glyph_source.glyph_source = g_src

        else:
            self._glyph = vtk_vectors()

        self._glyph.glyph.glyph.scale_factor = self._scale_factor
        #self._glyph.trait_modified = True
        self._glyph.glyph.glyph_source.glyph_position = self._position
        if not self._scaling:
            self._glyph.glyph.glyph.scaling = False
        self._glyph.actor.mapper.interpolate_scalars_before_mapping = True
        if self._glossy:
            self._glyph.actor.actor.property.specular = 0.7 * self._glossy
            self._glyph.actor.actor.property.specular_power = 11 * self._glossy
        if not self._shading:
            self._glyph.actor.actor.property.lighting = False

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

    def _update_lut_range(self, **kwargs):
        lut_range = np.array(self.data_range)
        self._module.vector_lut_manager.data_range = lut_range
        self._module.vector_lut_manager.use_default_range = False

    def _update_glyph(self, **kwargs):
        self._scale_factor = kwargs.get('scale_factor', self._scale_factor)
        self._glyph.glyph.glyph.scale_factor = self._scale_factor

    def _update_data(self, **kwargs):
        self._vectors = kwargs.get('vectors', self._vectors)
        self._scalars = kwargs.get('scalars', self._scalars)
        norm = np.linalg.norm(self._vectors, axis=1)
        mask = np.where(norm > 0)[0]
        self.data_range = [0, np.max(norm)]
        if self._geometry is not None:
            if self._geometry.type is 'Mesh':
                #self._geometry = self._anchor
                if self._anchor_mode == 'vertex':
                    anchor = self._geometry.vertices
                elif self._anchor_mode == 'face':
                    anchor = self._geometry.face_barycenters()
                elif self._anchor_mode == 'edge':
                    anchor = self._geometry.edge_mid_points()
                anchor = anchor[mask]
                self._data.points=anchor
            elif self._anchor.type is 'Points':
                self._geometry = self._anchor
                anchor = self._geometry.vertices
                anchor = anchor[mask]
                self._data.points=anchor
        elif 'anchor' in kwargs:
            anchor = kwargs['anchor']
            anchor = anchor[mask]
            self._data.points=anchor
        vectors = self._vectors[mask]
        if self._scalars is not None:
            self._scalars = np.array(self._scalars)
            scalars = self._scalars[mask]
            self._data.point_data.scalars = scalars
            self._data.point_data.scalars.name = 'scalars'
        self._data.point_data.vectors = vectors
        self._data.point_data.vectors.name = self.name
        self._data.modified()

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    def update(self, **kwargs):
        self._update_data(**kwargs)
        self._update_lut_range(**kwargs)
        self._update_glyph(**kwargs)
        self.source.update()

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------