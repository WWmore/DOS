""" A viewer for mlab scene. Adds a button to open up the engine.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2008, Enthought, Inc.
# License: BSD Style.

# Standard library imports
from os.path import join

import numpy as np

# Enthought library imports
from tvtk.tools.ivtk import IVTK
from tvtk.pyface.api import DecoratedScene
from traits.api import Callable
from pyface.api import ImageResource
from pyface.action.api import Action, Group
from pyface.resource.api import resource_path

# Local imports
from mayavi.core.common import error
from mayavi.preferences.api import set_scene_preferences, \
    get_scene_preferences


###############################################################################
# A decorated scene with an additional button.
###############################################################################
class GeolabScene(DecoratedScene):
    image_search_path = [join(resource_path(), 'images'), ]

    _position = None

    ##########################################################################
    # Non-public interface.
    ##########################################################################
    def show_engine(self):
        """ Open the engine view corresponding to the engine of the
            scene.
        """
        from mayavi.core.registry import registry
        from mayavi.core.ui.engine_rich_view import EngineRichView
        engine = None
        try:
            engine = registry.find_scene_engine(self)
        except TypeError:
            error('This scene is not managed by Mayavi')
        return EngineRichView(engine=engine).scene_editing_view(scene=self)

    ######################################################################
    # Trait handlers.
    ######################################################################
    def _actions_default(self):
        actions = [Group(
            Action(tooltip="View the pipeline",
                   image=ImageResource('img/new/mayavi.png'),
                   on_perform=self.show_engine,
                   ),
        ),
            Group(
                Action(tooltip="Save current view",
                       image=ImageResource('img/new/saveview.png'),
                       on_perform=self._save_view,
                       ),
                Action(tooltip="Apply saved view",
                       image=ImageResource('img/new/setview.png'),
                       on_perform=self._set_view,
                       ),
            ),
        ]
        actions.extend(DecoratedScene._actions_default(self))
        return actions

    def _save_view(self):
        self._position = self._get_position()
        print('position = ' + str(self._position))

    def _get_position(self):
        cc = self.camera
        p = [cc.position[0], cc.position[1], cc.position[2]]
        p.extend([cc.view_up[0], cc.view_up[1], cc.view_up[2]])
        p.extend([cc.focal_point[0], cc.focal_point[1], cc.focal_point[2]])
        p.extend([cc.view_angle])
        p.extend([cc.clipping_range[0], cc.clipping_range[1]])
        p.extend([self.parallel_projection])
        return p

    def _set_view(self, position=None):
        cc = self.camera
        p = self._position
        if position is not None:
            p = position
        if p is not None:
            cc.position = np.array([p[0], p[1], p[2]])
            cc.view_up = np.array([p[3], p[4], p[5]])
            cc.focal_point = np.array([p[6], p[7], p[8]])
            cc.view_angle = np.float(p[9])
            cc.clipping_range = np.array([p[10], p[11]])
            try:
                cc.parallel_projection = p[12]
            except:
                pass
            self.render()


def mayavi_scene_factory(parent):
    """A mayavi scene factory that creates a scene with preferences
    appropriately set."""
    p = get_scene_preferences()
    s = GeolabScene(parent, stereo=p['stereo'])
    set_scene_preferences(s, p)
    return s


###############################################################################
# A viewer making use of the MayaviScene
###############################################################################
class GeolabViewer(IVTK):
    """ A viewer window for mlab.
    """

    _scene_factory = Callable(mayavi_scene_factory)

    def _size_default(self):
        return 600, 300


def viewer_factory(size=(400, 350)):
    viewer = GeolabViewer()
    viewer.menu_bar_manager = None
    viewer.size = size
    viewer.open()
    return viewer


if __name__ == '__main__':
    from mayavi.tools.show import show

    viewer_factory()
    show()
