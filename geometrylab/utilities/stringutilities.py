# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import sys

import os

#------------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------


def print_progress(iteration, max_iteration):
    out = '|' + '#'*iteration + '.'*(max_iteration-iteration) + '|\r'
    sys.stdout.write(out)
    sys.stdout.flush()

def make_filepath(file_name, ext='txt', overwrite=False):
    folder = os.getcwd()
    name = os.path.join(folder, file_name)
    path = ('{}.{}').format(name, ext)
    if not overwrite:
        n = 1
        while os.path.exists(path):
            path = ('{}_({}).{}').format(name, n, ext)
            n += 1
    return path