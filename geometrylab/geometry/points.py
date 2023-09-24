#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from scipy import spatial

# -----------------------------------------------------------------------------

'''polyine.py: The Points data structure class'''

__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                                      Points
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Points(np.ndarray):

    def __new__(cls, input_array, info=None):
        obj = np.asarray(input_array).view(cls)
        obj.name = 'points'
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.name = getattr(obj, 'name', None)

    #--------------------------------------------------------------------------
    #                               Attributes
    #--------------------------------------------------------------------------

    @property
    def type(self):
        return 'Points'

    @property
    def V(self):
        return self.shape[0]

    @property
    def vertices(self):
        return self

    @property
    def points(self):
        return self


    #--------------------------------------------------------------------------
    #                                  CLOSENESS
    #--------------------------------------------------------------------------

    def closest_points(self, points, make_tree=False):
        KDtree = spatial.c_kDTree(self.vertices)
        closest = KDtree.query(points)[1]
        return closest

