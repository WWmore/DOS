#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from tvtk.api import tvtk

from mayavi.filters.poly_data_normals import PolyDataNormals

from mayavi.sources.vtk_data_source import VTKDataSource

from mayavi.modules.surface import Surface

from mayavi.core.api import ModuleManager

# -----------------------------------------------------------------------------

from geometrylab.vtkplot import plotutilities

from geometrylab.vtkplot import glyphs

# -----------------------------------------------------------------------------

'''edgesource.py: The mesh edges plot source class'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                   Edges
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class MeshVectors:

    def __init__(self,

                 vectors,

                 anchor,

                 position='head',

                 offset=0,

                 data=None,

                 lut_mode='blue-red',

                 reverse_lut=False,

                 lut_range=None,

                 color=None,

                 name='vectors',

                 scalars=None,

                 scale_factor=1,

                 resolution=20,

                 shaft_radius=None,

                 tip_length=None,

                 tip_radius=None,

                 smooth=True,

                 glossy=0,

                 specular=1):
        #----------------------------------------------------------------------
        try:
            anchor = anchor.vertices
        except AttributeError:
            pass
        mesh = glyphs.Vectors(vectors, anchor,shaft_radius, tip_length,
                              tip_radius, position, offset, resolution)
        #--------------
        self.name = name
        N = len(vectors)
        #------------------------ set tvtk_poly_data ----------------------------
        cell_types = np.array([tvtk.Triangle().cell_type,
                               tvtk.Quad().cell_type,
                               tvtk.Polygon().cell_type])
        cell_types = cell_types[mesh.types]
        cells = np.array([mesh.cells])
        cell_array = tvtk.CellArray()
        cell_array.set_cells(mesh.F, cells)
        ug = tvtk.UnstructuredGrid(points=mesh.vertices)
        ug.set_cells(cell_types, cell_array)
        source = ug
        if data is not None:
            scalars = scalars = np.repeat(data,mesh.F//N+1)
            self.lut_range = [np.min(data),np.max(data)]
        elif color is None:
            L = mesh.lengths
            scalars = np.repeat(L,mesh.F//N+1)
            self.lut_range = [np.min(L),np.max(L)]
        else:
            scalars = np.zeros((mesh.vertices.shape[0]))
        source.point_data.scalars = scalars
        source.point_data.scalars.name = name
        self.data = name
        #---------------------------- surface ---------------------------------
        surface = Surface()
        #---------------------------- smooth ----------------------------------
        if smooth:
            normals = PolyDataNormals()
            surface.actor.actor.property.specular = 0.1 #0.1
            surface.actor.actor.property.specular_power = 7 #7
        else:
            normals = None
        if glossy != 0:
            surface.actor.actor.property.specular = 0.7 * glossy
            surface.actor.actor.property.specular_power = 17 * specular
        #---------------------------- module ----------------------------------
        module = ModuleManager()
        if lut_range is not None:
            self.lut_range = lut_range
        if color is None:
            module.scalar_lut_manager.lut_mode = lut_mode
            module.scalar_lut_manager.use_default_range = False
            module.scalar_lut_manager.data_range = np.array(self.lut_range)
        else:
            lut = plotutilities.lut_table(color,1)
            module.vector_lut_manager.lut.table = lut
            self.lut_range = None
        module.add_child(surface)
        #----------------------------------------------------------------------
        if normals is None:
            src = VTKDataSource(data = source)
            src.add_child(module)
        else:
            src = VTKDataSource(data = source)
            normals.add_child(module)
            src.add_child(normals)
        self.source = src