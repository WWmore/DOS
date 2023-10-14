#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import os, sys

from traits.api import Button, Bool, Int, on_trait_change

from traitsui.api import View, Item, VGroup

import numpy as np

#------------------------------------------------------------------------------

path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# print(path)
sys.path.append(path)
                
import geometrylab as geo

#------------------------------------------------------------------------------

''' Sample component for GeoLab.
__author__ = 'Davide Pellis'

Hui note: archgeolab refers to this class; if want to make new one, refer here!
'''




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                             SAMPLE COMPONENT CLASS
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


class Paneling(geo.gui.GeolabComponent):

    name = 'Paneling'

    #--------------------------------------------------------------------------
    #                                 Traits
    #--------------------------------------------------------------------------

    set_something_input = Int(10, label = 'number of points')

    do_something_button = Button(label = 'do something')

    move_and_do_something_check = Bool(False, label = 'move and do something')

    generate_points_button = Button(label = 'generate points')

    generate_vectors_button = Button(label = 'generate vectors')

    make_dual_mesh_button = Button(label = 'make dual mesh')

    #--------------------------------------------------------------------------
    #                              Component View
    #--------------------------------------------------------------------------

    view = View(
                VGroup(
                       Item('set_something_input'),
                       Item('move_and_do_something_check'),
                       Item('do_something_button',
                            show_label=False),
                       Item('generate_points_button',
                            show_label=False),
                       Item('generate_vectors_button',
                            show_label=False),
                       Item('make_dual_mesh_button',
                            show_label=False),
                       show_border=True,
                       ),
                resizable=True,
                )

    #--------------------------------------------------------------------------
    #                                Attributes
    #--------------------------------------------------------------------------

    counter = 0

    #--------------------------------------------------------------------------
    #                            Special Methods
    #--------------------------------------------------------------------------
    '''These are overwriting of special methods, fired by geolab'''

    def geolab_settings(self):
        self.geolab.height = 800
        self.geolab.width = 900
        self.geolab.side_width = 300

    def object_open(self, file_name, geometry):
        name = ('mesh_{}').format(self.counter)
        self.geolab.add_object(geometry, name=name)
        self.counter += 1
        obj = self.geolab.last_object
        if obj.geometry.type == 'Mesh':
            self.initialize_mesh_plot_functions(obj)

    def object_change(self):
        '''Function launched before an object selection change'''
        pass

    def object_changed(self):
        '''Function launched after an object selection change or after a
        topology change of the current object'''
        pass

    def object_save(self, file_name):
        '''Function launched after a mesh is saved. The file_name is generated
        by the main application'''
        pass

    def set_state(self, state):
        if state != 'move_and_do_something':
            self.move_and_do_something_check = False

    #--------------------------------------------------------------------------
    #                              Plot Functions
    #--------------------------------------------------------------------------

    def initialize_mesh_plot_functions(self, obj):

        def plot_something_on_faces():
            data = np.random.random(obj.geometry.F)
            obj.plot_faces(face_data = data,
                           color = 'blue-red',
                           lut_range = '-:+',
                           smooth = True)
        obj.add_face_callback(plot_something_on_faces,
                              name='something on faces')

        obj.face_plot = 'something on faces'

        def plot_something_on_vertices():
            data = obj.geometry.vertices[:,2]
            obj.plot_faces(vertex_data = data,
                           color = 'blue-red',
                           lut_range = '-:+',
                           smooth = True)
        obj.add_face_callback(plot_something_on_vertices,
                              name='something on vertices')

        def plot_something_on_edges():
            data = np.random.random(obj.geometry.E)
            obj.plot_edges(edge_data = data,
                           color = 'blue-red',
                           lut_range = '-:+')
        obj.add_edge_callback(plot_something_on_edges,
                            name='something on edges')

        obj.edge_plot = 'something on edges'

    def initialize_points_plot_functions(self, obj):

        def plot_something_on_points():
            data = np.random.random(obj.geometry.V)
            obj.plot_vertices(vertex_data = data,
                              color = 'blue-red',
                              lut_range = '-:+')
        obj.add_vertices_callback(plot_something_on_points,
                                  name='something on vertices')

        obj.vertices_plot = 'something on vertices'

    #--------------------------------------------------------------------------
    #                              Traits Change
    #--------------------------------------------------------------------------

    @on_trait_change('set_something_input')
    def set_something(self):
        '''The action fired when set_something is changed'''
        self.generate_points()

    @on_trait_change('do_something_button')
    def do_something(self):
        '''The action fired when do-something is pressed'''
        pass

    @on_trait_change('generate_points_button')
    def generate_points(self):
        '''Make a points geometry'''
        n = self.set_something_input
        P = geo.geometry.Points(np.random.random((n, 3))*10 + 4)

        '''Add the points to geolab. Since the name is not changing,
        the object is replaced'''
        self.geolab.add_object(P, name='points')

        obj = self.geolab.get_object('points')

        self.initialize_points_plot_functions(obj)

    @on_trait_change('generate_vectors_button')
    def generate_vectors(self):
        '''Generates random vectors on the current object vertices'''
        obj = self.geolab.current_object

        '''Plot the vectors over the object vertices, adding a callback.
        This links the vectors to the object in geometry updates'''
        def plot_vectors():
            vec = np.random.random((obj.geometry.V, 3))*20
            obj.plot_vectors(vectors = vec,
                             color = 'blue-red',
                             name = 'vectors')
        obj.add_plot_callback(plot_vectors, name = 'vectors')
        obj.update_plot()


    @on_trait_change('make_dual_mesh_button')
    def make_dual_mesh(self):
        '''Get the current object and check if it is a mesh'''
        obj = self.geolab.current_object
        if obj.geometry.type != 'Mesh':
            return

        '''Make the dual geometry'''
        D = obj.geometry.dual_mesh()
        D.move([0,12,0])
        G = geo.optimization.Gridshell()
        G.import_mesh(D)

        '''Give a name to the dual, unique for each primal mesh, and
        add it to geolab and initialize its plot. In this way, subsequent
        dual meshes of the same primal will be replaced.'''
        name = obj.name + '_dual'
        self.geolab.add_object(G, name = name)
        self.geolab.add_object(G, name = name + '_copy', scene = 'scene_1') ##Hui: to debug
        self.initialize_mesh_plot_functions(self.geolab.get_object(name))
        self.initialize_mesh_plot_functions(self.geolab.get_object(name+'_copy'))
        #self.geolab.get_scene('scene_1').parallel_projection()
        #self.geolab.get_scene('scene_1').z_view()
        self.geolab.get_scene('scene_1').interaction_2d()


    @on_trait_change('move_and_do_something_check')
    def move_and_do_something(self):
        '''This set the state avoiding conflicts'''
        self.handler.set_state('move_and_do_something')

        '''Get the selected object'''
        obj = self.geolab.current_object

        if self.move_and_do_something_check:

            '''The callback for point interaction, as argument is passed the
            selected vertex index'''
            def vertex_callback():

                '''The index of the moving vertex'''
                selected = obj.selected_vertices

                '''The coordinates of the point'''
                point = obj.geometry.vertices[selected]

                '''Plot a new point as a source (not selectable)'''
                obj.plot_points(point/2,
                                radius=0.3,
                                color='r',
                                name='half_point')

                '''Update the plot of the object and its sources'''
                #obj.update_plot()
                self.geolab.update_plot()

            '''Start the moving vertex loop'''
            obj.move_vertices(vertex_callback)

        else:

            '''Stop the moving vertex loop'''
            obj.move_vertices_off()

            '''Hide the points'''
            obj.hide('half_point')

            obj.update_plot()


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#                               COMPONENT TEST
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#                     File selection from the file folder
#------------------------------------------------------------------------------

path = os.path.dirname(os.path.abspath(__file__))

file_name = path + '/tri_dome.obj'

#------------------------------------------------------------------------------
#                     Main code to run the application
#------------------------------------------------------------------------------

if __name__ == '__main__':

    '''Instantiate the sample component'''
    component = Paneling()

    '''Instantiate the main geolab application'''
    GUI = geo.gui.GeolabGUI()

    '''Add the component to geolab'''
    GUI.add_component(component)

    '''Open an obj file'''
    GUI.open_obj_file(file_name)

    '''Start geolab main loop'''
    GUI.start()


