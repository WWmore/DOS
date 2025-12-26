#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from traits.api import Instance, String,on_trait_change, Bool, Int, List

from traitsui.api import View, Item, HSplit, ListEditor, Action, ToolBar, Separator, Controller

from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel

from pyface.image_resource import ImageResource

#------------------------------------------------------------------------------

from geometrylab.gui.scenemanager import SceneManager

from geometrylab.gui.multiscenemanager import MultiSceneManager

from geometrylab.gui.geolabscene import GeolabScene

from geometrylab.gui.geolabcomponent import GeolabComponent

from geometrylab.gui.handler import Handler

from geometrylab.optimization.gridshell import Gridshell

from geometrylab.gui.tools import IncrementalRemesh, CornerTolerance,\
                                  Loads, SaveMesh



# -----------------------------------------------------------------------------

__author__ = 'Davide Pellis'


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                                     Handler
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class GlHandler(Controller):

    def close(self, info, is_ok):
        info.object._closed = True
        Controller.close(self ,info, is_ok)
        return True

    #--------------------------------------------------------------------------
    #                                  Tools
    #--------------------------------------------------------------------------

    def open_plot_manager(self):
        self.info.object.open_plot_manager()

    def open_corner_tolerance(self):
        self.info.object.open_corner_tolerance()

    def open_save_mesh(self):
        self.info.object.open_save_mesh()

    def open_remesh(self):
        self.info.object.open_remesh()

    def open_loads(self):
        self.info.object.open_loads()


    #--------------------------------------------------------------------------
    #                                 Selection
    #--------------------------------------------------------------------------

    def background_switch(self):
        self.info.object.background_switch()

    def select_object(self):
        self.info.object.select_object()

    def select_vertices(self):
        self.info.object.select_vertices()

    def select_edges(self):
        self.info.object.select_edges()

    def select_faces(self):
        self.info.object.select_faces()

    def select_boundary_vertices(self):
        self.info.object.select_boundary_vertices()

    def move_vertices(self):
        self.info.object.move_vertices()

    #--------------------------------------------------------------------------
    #                                 Mesh state
    #--------------------------------------------------------------------------

    def reset_mesh(self):
        self.info.object.reset_mesh()

    def set_reference(self):
        self.info.object.set_reference()

    #--------------------------------------------------------------------------
    #                                 Mesh Edit
    #--------------------------------------------------------------------------

    def flip_edges(self):
        self.info.object.flip_edges()

    def collapse_edges(self):
        self.info.object.collapse_edges()

    def split_edges(self):
        self.info.object.split_edges()

    def catmull_clark(self):
        self.info.object.catmull_clark()

    def loop(self):
        self.info.object.loop()

    #--------------------------------------------------------------------------
    #                              Set conditions
    #--------------------------------------------------------------------------

    def fix_vertices(self):
        self.info.object.fix_vertices()

    def unfix_vertices(self):
        self.info.object.unfix_vertices()

    def constrain_vertices(self):
        self.info.object.constrain_vertices()


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                       The Base Graphic User Interface
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


class GeolabGUI(MultiSceneManager):

    _handler = GlHandler()

    _closed = Bool(True)

    _current_object = String('none')

    _components = []

    scene_model_0 = Instance(MlabSceneModel, ())

    scene_model_1 = Instance(MlabSceneModel, ())

    scene_model_2 = Instance(MlabSceneModel, ())

    scene_model_3 = Instance(MlabSceneModel, ())

    components = List(GeolabComponent)

    remesh_tool = IncrementalRemesh()

    corner_tolerance_tool = CornerTolerance()

    loads_tool = Loads()

    save_mesh_tool = SaveMesh()

    _background_switch_counter = Int(1)

    background_switch_button = Action(action = 'background_switch',
                               image = ImageResource('img/new/background.png'),
                               style = 'push',
                               tooltip = 'Switch background color',
                               show_label = False)

    save_mesh_button = Action(action = 'open_save_mesh',
                       image = ImageResource('img/new/savemesh.png'),
                       style = 'push',
                       tooltip = 'Save mesh',
                       show_label = False)

    plot_manager_button = Action(action = 'open_plot_manager',
                          image=ImageResource('img/new2/plotmanager.png'),
                          style = 'push',
                          tooltip = 'Object plot settings',
                          enabled_when = ('_current_object == \
                                          "Mesh_plot_manager" \
                                          or _current_object == \
                                          "Points_plot_manager"'),
                          show_label=False)

    corner_tolerance_button = Action(action = 'open_corner_tolerance',
                              image = ImageResource('img/new2/corners.png'),
                              style = 'push',
                              tooltip = 'Set corner tolerance',
                              enabled_when = ('_current_object == \
                                              "Mesh_plot_manager"'),
                              show_label=False)

    select_object_button = Action(action = 'select_object',
                           image = ImageResource('img/new2/selectobject.png'),
                           style = 'toggle',
                           tooltip = 'Select object',
                           show_label = False)

    flip_edges_button = Action(action = 'flip_edges',
                        image = ImageResource('img/new2/flip.png'),
                        style = 'toggle',
                        tooltip = 'Flip mesh edges',
                        enabled_when = ('_current_object == \
                                        "Mesh_plot_manager"'),
                        show_label=False)

    split_edges_button = Action(action='split_edges',
                         image = ImageResource('img/new2/split.png'),
                         style = 'toggle',
                         enabled_when = ('_current_object == \
                                         "Mesh_plot_manager"'),
                         tooltip = 'Split mesh edges',
                         show_label = False)

    collapse_edges_button = Action(action='collapse_edges',
                            image=ImageResource('img/new2/collapse.png'),
                            style='toggle',
                            enabled_when =('_current_object == \
                                           "Mesh_plot_manager"'),
                            tooltip = 'Collapse mesh edges',
                            show_label=False)

    catmull_clark_button = Action(action = 'catmull_clark',
                           image = ImageResource('img/new2/catmullclark.png'),
                           style = 'push',
                           enabled_when = ('_current_object == \
                                           "Mesh_plot_manager"'),
                           tooltip = 'Catmull-Clark subdivision',
                           show_label = False)

    loop_button = Action(action = 'loop',
                  image = ImageResource('img/new2/loop.png'),
                  style = 'push',
                  enabled_when = ('_current_object == "Mesh_plot_manager"'),
                  tooltip = 'Loop subdivision',
                  show_label = False)

    remesh_button = Action(action = 'open_remesh',
                    image = ImageResource('img/new2/remesh.png'),
                    style = 'push',
                    enabled_when = ('_current_object == "Mesh_plot_manager"'),
                    tooltip = 'Incremental remesh',
                    show_label = False)

    move_vertices_button = Action(action = 'move_vertices',
                           image = ImageResource('img/new2/movevertex.png'),
                           style = 'toggle',
                           tooltip = 'Move vertices',
                           enabled_when = ('_current_object == \
                                           "Mesh_plot_manager" \
                                           or _current_object == \
                                           "Points_plot_manager"'),
                           show_label=False)

    select_vertices_button = Action(action='select_vertices',
                             image=ImageResource('img/new2/selectvertices.png'),
                             style='toggle',
                             tooltip = 'Select vertices',
                             enabled_when = ('_current_object == \
                                             "Mesh_plot_manager" \
                                             or _current_object == \
                                             "Points_plot_manager"'),
                             show_label=False)

    select_edges_button = Action(action='select_edges',
                          image=ImageResource('img/new2/selectedges.png'),
                          style='toggle',
                          enabled_when =('_current_object=="Mesh_plot_manager"'),
                          tooltip = 'Select edges',
                          show_label=False)

    select_faces_button = Action(action='select_faces',
                          image=ImageResource('img/new2/selectfaces.png'),
                          style='toggle',
                          enabled_when =('_current_object == \
                                         "Mesh_plot_manager"'),
                          tooltip = 'Select faces',
                          show_label=False)

    select_boundary_vertices_button = Action(action='select_boundary_vertices',
                                   image=ImageResource('img/new2/boundary.png'),
                                   style='toggle',
                                   enabled_when =('_current_object == \
                                                  "Mesh_plot_manager"'),
                                   tooltip = 'Select boundary vertices',
                                   show_label=False)

    fix_vertices_button = Action(action='fix_vertices',
                          image=ImageResource('img/new/fixvertices.png'),
                          style='push',
                          tooltip = 'Fix selected vertices',
                          show_label=False)

    unfix_vertices_button = Action(action='unfix_vertices',
                            image=ImageResource('img/new/unfixvertices.png'),
                            style='push',
                            tooltip = 'Unfix selected vertices',
                            show_label=False)

    constrain_vertices_button = Action(action='constrain_vertices',
                                image=ImageResource('img/new/constrain.png'),
                                style='push',
                                tooltip = 'Constrain selected vertices',
                                show_label=False)

    unconstrain_vertices_button = Action(action='constrain_vertices',
                                image=ImageResource('img/new/unconstrain.png'),
                                style='push',
                                tooltip = 'Release selected vertices',
                                show_label=False)

    loads_button = Action(action='open_loads',
                   image=ImageResource('img/new/applyforce.png'),
                   style='push',
                   tooltip = 'loads_tool',
                   show_label=False)

    reset_mesh_button = Action(action='reset_mesh',
                        image=ImageResource('img/new/resetmesh.png'),
                        style='push',
                        tooltip = 'Reset mesh',
                        show_label=False)

    set_reference_button = Action(action='set_reference',
                           image=ImageResource('img/new/setreference.png'),
                           style='push',
                           tooltip = 'Set as reference mesh',
                           show_label=False)

    __toolbar = [Separator(),
                 plot_manager_button,

                 Separator(),
                 select_object_button,
                 select_vertices_button,
                 select_edges_button,
                 select_faces_button,
                 move_vertices_button,
                 select_boundary_vertices_button,
                 corner_tolerance_button,

                 Separator(),
                 fix_vertices_button,
                 unfix_vertices_button,
                 constrain_vertices_button,
                 unconstrain_vertices_button,
                 loads_button,

                 Separator(),
                 collapse_edges_button,
                 flip_edges_button,
                 split_edges_button,
                 catmull_clark_button,
                 loop_button,
                 remesh_button,

                 Separator(),
                 save_mesh_button,
                 reset_mesh_button,
                 set_reference_button,

                 Separator(),
                 background_switch_button,
                 ]

    toolbar = { 'handler': _handler,
                'resizable' : True,
                'title' : 'GeometryLab',
                'icon' : ImageResource('img/new2/logo3.png'),
                'toolbar' : ToolBar(*__toolbar,
                                    show_labels=False,
                                    image_size=(16,16)),
              }



    tabs = [Item('components',
                 editor=ListEditor(use_notebook=True, page_name='.name'),
                 style ='custom',
                 width = 1,
                 resizable = False,
                 show_label=False)
           ]

    __view3D = []

    __windows = []

    height = Int(300)

    width = Int(500)

    side_width = Int(400)

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    def __init__(self):
        MultiSceneManager.__init__(self)

        self.handler = Handler()

        self.__scenes = ['scene_0']

        self.__geometries = []

        self.__object_open_callbacks = []

        self.__initialize_plot = True

        self._scene_model = self.scene_model_0

        self.handler.add_state_callback(self.set_state)

    @property
    def geometries(self):
        return self.__geometries

    #--------------------------------------------------------------------------
    #                             Building up
    #--------------------------------------------------------------------------

    def start(self):
        self.make_scenes()
        if len(self._components) > 0 and len(self.__windows) > 0:
            self.components = self._components
            view = View(HSplit(self.tabs,
                               self.__view3D,
                               self.__windows),
                        **self.toolbar)
        elif len(self._components) > 0:
            self.components = self._components
            view = View(HSplit(self.tabs,
                               self.__view3D),
                        **self.toolbar)
        else:
            view = View(self.__view3D, **self.toolbar)
        self._closed = False
        for component in self._components:
            component.initialize_plot()
        if self.__initialize_plot:
            for geometry in self.__geometries:
                self.add(geometry)
            for key in self._objects:
                self._objects[key].update_plot()
        for args in self.__object_open_callbacks:
            self.object_open(*args)
        self.object_changed()
        self.update_plot()
        self.configure_traits(view=view)

    def add_component(self, component):
        component.geolab = self
        self._components.append(component)
        component.geolab_settings()
        self.__initialize_plot = False

    def add_scene(self, name):
        if name not in self.__scenes:
            self.__scenes.append(name)

    def make_scenes(self):
        index = 0
        for key in self.__scenes:
            if index == 0:
                scene = self.scene_model_0
                self.__view3D = [Item('scene_model_0',
                                 editor = SceneEditor(scene_class=GeolabScene),
                                 show_label = False,
                                 resizable = True,
                                 height = self.height,
                                 width = self.width)]
            else:
                if index == 1:
                    scene = self.scene_model_1
                elif index == 2:
                    scene = self.scene_model_2
                elif index == 3:
                    scene = self.scene_model_3
                else:
                    return
                name = ('scene_model_{}').format(index)
                self.__windows.append(Item(name,
                                 editor = SceneEditor(scene_class=GeolabScene),
                                 show_label = False,
                                 resizable = True,
                                 height = self.height,
                                 width = self.side_width))
            self.add_scene_model(scene, key)
            index += 1

    #--------------------------------------------------------------------------
    #                             Standard Methods
    #--------------------------------------------------------------------------

    def set_state(self, name):
        self.select_off()
        if name != 'select_object':
            self.select_object_button.checked = False
        if name != 'flip_edges':
            self.flip_edges_button.checked = False
        if name != 'split_edges':
            self.split_edges_button.checked = False
        if name != 'collapse_edges':
            self.collapse_edges_button.checked = False
        if name != 'move_vertices':
            self.move_vertices_button.checked = False
        if name != 'select_vertices':
            self.select_vertices_button.checked = False
        if name != 'select_faces':
            self.select_faces_button.checked = False
        if name != 'select_edges':
            self.select_edges_button.checked = False
        if name != 'select_boundary_vertices':
            self.select_boundary_vertices_button.checked = False

    def object_changed(self):
        SceneManager.object_changed(self)
        self.close_tools()
        self._current_object = self.current_object_type

    #--------------------------------------------------------------------------
    #                                Background
    #--------------------------------------------------------------------------

    def background_switch(self):
        if self._background_switch_counter == 1:
            self.set_background((1, 1, 1))
            self._background_switch_counter = 0
        else:
            self.set_background((0.5, 0.5, 0.5))
            self._background_switch_counter = 1

    #--------------------------------------------------------------------------
    #                                  Open
    #--------------------------------------------------------------------------

    def open_obj_file(self, file_name):
        mesh = Gridshell()
        mesh.read_obj_file(file_name)
        self.__geometries.append(mesh)
        if not self._closed:
            self.object_open(file_name, mesh)
        else:
            self.__object_open_callbacks.append((file_name, mesh))

    def open_geometry(self, geometry):
        self.__object_open_callbacks.append(('M', geometry))

    def reset_mesh(self):
        self.handler.set_state(None)
        self.current_object.geometry.reset()
        self.current_object.update_plot()
        self.object_changed()

    def set_reference(self):
        self.current_object.geometry.set_reference()

    #--------------------------------------------------------------------------
    #                                  Tools
    #--------------------------------------------------------------------------

    def open_plot_manager(self):
        self.selection_off()
        self.current_object.start()

    def open_corner_tolerance(self):
        self.handler.set_state(None)
        self.corner_tolerance_tool.scenemanager = self
        self.corner_tolerance_tool.start()

    def open_save_mesh(self):
        self.save_mesh_tool.scenemanager = self
        self.save_mesh_tool.start()

    def open_remesh(self):
        self.handler.set_state(None)
        self.remesh_tool.scenemanager = self
        self.remesh_tool.start()

    def open_loads(self):
        self.loads_tool.scenemanager = self
        self.loads_tool.start()

    @on_trait_change('_closed')
    def close_tools(self):
        try:
            self.current_object.close()
        except:
            pass
        self.corner_tolerance_tool.close()
        self.remesh_tool.close()
        self.loads_tool.close()
        self.save_mesh_tool.close()

    #--------------------------------------------------------------------------
    #                                Selection
    #--------------------------------------------------------------------------

    def selection_off(self):
        MultiSceneManager.select_off(self)
        self.select_boundary_vertices_button.checked = False
        self.select_edges_button.checked = False
        self.select_object_button.checked = False
        self.select_vertices_button.checked = False
        self.move_vertices_button.checked = False

    def select_object(self):
        if self.select_object_button.checked:
            self.handler.set_state('select_object')
            super(GeolabGUI, self).select_object()
        else:
            self.select_off()

    def select_vertices(self):
        if self.select_vertices_button.checked:
            self.handler.set_state('select_vertices')
            self.current_object.select_vertices()
        else:
            self.current_object.select_vertices_off()

    def select_edges(self):
        if self.select_edges_button.checked:
            self.handler.set_state('select_edges')
            self.current_object.select_edges()
        else:
            self.current_object.select_edges_off()

    def select_faces(self):
        if self.select_faces_button.checked:
            self.handler.set_state('select_faces')
            self.current_object.select_faces()
        else:
            self.current_object.select_faces_off()

    def select_boundary_vertices(self):
        if self.select_boundary_vertices_button.checked:
            self.handler.set_state('select_boundary_vertices')
            self.current_object.select_boundary_vertices()
        else:
            self.current_object.select_boundary_vertices_off()

    def move_vertices(self):
        if self.move_vertices_button.checked:
            self.handler.set_state('move_vertices')
            def callback():
                self.current_object.update_plot()
            self.current_object.move_vertices(callback)
        else:
            self.current_object.move_vertices_off()

    #--------------------------------------------------------------------------
    #                                Mesh Edit
    #--------------------------------------------------------------------------

    def flip(self, edge_index):
        self.current_object.mesh.flip_edge(edge_index)
        self.object_changed()
        self.current_object.update_plot()

    def collapse(self, edge_index):
        self.current_object.mesh.collapse_edge(edge_index)
        self.object_changed()
        self.current_object.update_plot()

    def split(self, edge_index):
        self.current_object.mesh.split_edge(edge_index)
        self.object_changed()
        self.current_object.update_plot()

    def flip_edges(self):
        if self.flip_edges_button.checked:
            self.handler.set_state('flip_edges')
            self.current_object.on_edge_selection(self.flip)

        else:
            self.current_object.select_edges_off()

    def collapse_edges(self):
        if self.collapse_edges_button.checked:
            self.handler.set_state('collapse_edges')
            self.current_object.on_edge_selection(self.collapse)
            self.object_changed()
        else:
            self.current_object.select_edges_off()

    def split_edges(self):
        if self.split_edges_button.checked:
            self.handler.set_state('split_edges')
            self.current_object.on_edge_selection(self.split)
        else:
            self.current_object.select_edges_off()

    def catmull_clark(self):
        self.handler.set_state(None)
        self.fix_edges_bug()
        self.current_object.mesh.catmull_clark()
        self.object_changed()
        self.current_object.update_plot()

    def loop(self):
        self.handler.set_state(None)
        self.fix_edges_bug()
        self.current_object.mesh.loop()
        self.object_changed()
        self.current_object.update_plot()

    def fix_edges_bug(self):
        obj = self.current_object
        if True:#obj.mesh.E > 2000:
            if obj.edge_plot != 'none':
                if obj.edge_plot != 'wireframe':
                    obj.remove_edges()

    #--------------------------------------------------------------------------
    #                              Set conditions
    #--------------------------------------------------------------------------

    def fix_vertices(self):
        selected = self.current_object.selected_vertices
        self.current_object.mesh.fix(selected)
        self.current_object.update_plot()

    def unfix_vertices(self):
        selected = self.current_object.selected_vertices
        self.current_object.mesh.unfix(selected)
        self.current_object.update_plot()

    def constrain_vertices(self):
        selected = self.current_object.selected_vertices
        self.current_object.mesh.constrain(selected)
        self.current_object.update_plot()


# -----------------------------------------------------------------------------

def view(objects, position=None):
    viewer = GeolabGUI()
    #viewer.position = position
    viewer.add(objects)
    viewer.start()
