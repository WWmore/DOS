# -*- coding: utf-8 -*-


#------------------------------------------------------------------------------
#                                   HANDLER
#------------------------------------------------------------------------------


class Handler(object):

    def __init__(self):

        self.handled = {}

        self.state_callbacks = []

        self.state_variables = {}


    def add_state_callback(self, callback):
        self.state_callbacks.append(callback)

    def set_state(self, state=None):
        for name in self.handled:
            try:
                self.handled[name]._set_state(state)
            except:
                pass
        for callback in self.state_callbacks:
            try:
                callback(state)
            except:
                pass

    def get_variable(self, name, default=None):
        try:
            variable = self.state_variables[name]
            return variable()
        except:
            variable = self.state_variables[name]
            return variable()

    def add_variable(self, name, value):
        self.state_variables[name] = value

