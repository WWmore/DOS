

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

#------------------------------------------------------------------------------

from traits.api import HasTraits

#------------------------------------------------------------------------------

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                  Component
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class GeolabComponent(HasTraits):

    name = 'component'

    def __init__(self):
        HasTraits.__init__(self)
        self.__geolab = None

    @property
    def geolab(self):
        return self.__geolab

    @geolab.setter
    def geolab(self, geolab):
        self.__geolab = geolab
        self.__geolab.add_callback('object_change', self.object_change)
        self.__geolab.add_callback('object_changed', self.object_changed)
        self.__geolab.add_callback('object_save', self.object_save)
        self.__geolab.add_callback('object_open', self.object_open)
        self.__geolab.handler.add_state_callback(self.set_state)

    @property
    def handler(self):
        return self.geolab.handler

    def geolab_settings(self):
        pass

    def initialize_plot(self):
        pass

    def set_state(self):
        pass

    def object_change(self):
        pass

    def object_changed(self):
        pass

    def object_save(self, file_name):
        pass

    def object_open(self, file_name, geometry):
        pass





