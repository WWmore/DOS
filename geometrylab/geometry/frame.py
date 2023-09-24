#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

from geometrylab import utilities

# -----------------------------------------------------------------------------

'''frame.py: The Frame data structure class'''

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                   Frame
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class Frame(object):

    def __init__(self, origin=[0,0,0], e1=[1,0,0], e2=[0,1,0], e3=[0,0,1]):

        self.origin = origin

        self.e1 = e1

        self.e2 = e2

        self.e3 = e3

    @property
    def N(self):
        return self._origin.shape[0]

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        origin = np.array(origin, dtype=np.float)
        if len(origin.shape) == 1:
            origin = np.array([origin])
        if origin.shape[1] != 3:
            raise ValueError('the origin of a frame must be a (n,3) vector')
        self._origin = origin

    @property
    def e1(self):
        return self._x

    @e1.setter
    def e1(self, e1):
        e1 = np.array(e1, dtype=np.float)
        if len(e1.shape) == 1:
            e1 = np.array([e1])
        if e1.shape[1] != 3:
            raise ValueError('the e1-axis of a frame must be a (n,3) vector')
        if e1.shape[0] != self.N and e1.shape[0] == 1:
            e1 = np.tile(e1, (self.N, 1))
        elif e1.shape[0] != self.N:
            raise ValueError('dimesnsion mismatch')
        self._x = e1

    @property
    def e2(self):
        return self._y

    @e2.setter
    def e2(self, e2):
        e2 = np.array(e2, dtype=np.float)
        if len(e2.shape) == 1:
            e2 = np.array([e2])
        if e2.shape[1] != 3:
            raise ValueError('the e2-axis of a frame must be a (n,3) vector')
        if e2.shape[0] != self.N and e2.shape[0] == 1:
            e2 = np.tile(e2, (self.N, 1))
        elif e2.shape[0] != self.N:
            raise ValueError('dimesnsion mismatch')
        self._y = e2

    @property
    def e3(self):
        return self._z

    @e3.setter
    def e3(self, e3):
        e3 = np.array(e3, dtype=np.float)
        if len(e3.shape) == 1:
            e3 = np.array([e3])
        if e3.shape[1] != 3:
            raise ValueError('the e3-axis of a frame must be a (n,3) vector')
        if e3.shape[0] != self.N and e3.shape[0] == 1:
            e3 = np.tile(e3, (self.N, 1))
        elif e3.shape[0] != self.N:
            raise ValueError('dimesnsion mismatch')
        self._z = e3
