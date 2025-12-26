#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from __future__ import unicode_literals

import numpy as np

from mayavi.core.api import ModuleManager

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def color_library():
    library = {

    'r' :              (255, 0, 0),

    'g' :              (0, 255, 0),

    'b' :              (0, 0, 255),

    'c' :              (0, 255, 255),

    'm' :              (255, 0, 255),

    'y' :              (255, 255, 0),

    'k' :              (0, 0, 0),

    'w' :              (255, 255, 255),

    'white' :          (255, 255, 255),

    'black' :          (0, 0, 0),

    'gray_40' :        (153, 153, 153),

    'gray_50' :        (128, 128, 128),

    'gray_60' :        (102, 102, 102),

    'magenta' :        (255, 0, 255),

    'cyan' :           (0, 255, 255),

    'azure' :          (0, 136, 255),

    'red_orange' :     (255, 120, 64),

    'orange' :         (255, 128, 0),

    'green' :          (0, 255, 128),

    'yellow' :         (255, 255 ,0),

    'golden_yellow' :  (255, 221, 0),

    'metallic_gold' :  (212, 175 , 55),

    'gold' :           (218, 165, 32),

    'aqua' :           (0, 255, 255),

    'cadet_blue' :     (95, 158, 160),

    'sandy_brown' :    (244, 164, 96),

    'white_blue' :     (204, 229, 255),

    'cornflower' :     (122, 163, 230),

    'chrome_yellow':   (255, 170, 0),

    'safety_orange':   (255, 119, 0),

    'carmine_red':     (255, 0, 51),

    'violet':          (153, 17, 170),

    'vivid_lime_green':(170, 221, 34),

    'smashed_pumpkin': (251, 109, 55),

    'icterine':        (251, 247, 94),

    'tiffany_blue':    (1, 198, 178),

    'sea_blue':        (0, 108, 150),

    'prussian_blue':   (0, 57, 89),

    'capri':           (1, 190, 254),

    'amber':           (255, 125, 0),

    'vivid_raspberry': (255, 0, 109),

    'spring_bud':      (173, 255, 2),

    'electric_violet': (143, 0, 255),

    'cyber':           (255, 211, 0),

    }
    return library

def palette_library():
    library = {

    'raimbow_1':       ['w', 'azure', 'safety_orange', 'violet',
                        'vivid_lime_green', 'carmine_red', (255, 186, 0),
                        (0, 60, 95)],

    'one_more_raimbow':['w', 'azure', 'safety_orange', 'violet',
                        'vivid_lime_green', 'carmine_red', 'chrome_yellow'],

    'lose_me_away':    ['w', 'smashed_pumpkin', 'icterine', 'tiffany_blue',
                        'sea_blue', 'prussian_blue'],

    'happy_colors':    ['w', 'capri', 'golden_yellow', 'amber',
                        'vivid_raspberry', 'spring_bud', 'electric_violet'],

    'milti_copat':     ['w', (30, 56, 136), (245, 230, 99), (156, 56, 72),
                        (71, 168, 189), (255, 173, 105)],

    'overtly_cautions':['w',  (26, 91, 109), (216, 201, 155), (216, 151, 60),
                        (189, 99, 47),  (39, 62, 71), (140, 179, 201),
                        (152, 150, 118)]

    }
    return library

def palette(key):
    try:
        palettes = palette_library()
        return palettes[key]
    except:
        return key

def map_color(key):
    library = color_library()
    try:
        color = library[key]
    except KeyError:
        if type(key) == tuple and len(key) == 3:
            color = key
        elif type(key) == str:
            color = (128 ,128, 128)
        elif key is None:
            color = None
        else:
            color = (128 ,128, 128)
    return color

def rgb_float_color(key):
    c = map_color(key)
    return (c[0]/255, c[1]/255, c[2]/255)

def lut_table(colors, opacity):
    '''
    elif isinstance(colors, unicode):
        colors = [colors]
    '''
    colors = palette(colors)
    if isinstance(colors, str):
        colors = [colors]
    elif isinstance(colors, tuple):
        colors = [colors]
    if isinstance(opacity, float):
        opacity = [opacity]
    if isinstance(opacity, int):
        opacity = [opacity]
    if len(colors) != len(opacity):
        opacity = np.repeat(opacity[0],len(colors))
    rgb = []
    for key in colors:
        rgb.append(map_color(key))
    lut = np.ones((256, 4)) * 255
    color = rgb[-1]
    alpha = opacity[-1]
    lut[:, 0] = color[0] * np.ones([256])
    lut[:, 1] = color[1] * np.ones([256])
    lut[:, 2] = color[2] *  np.ones([256])
    lut[:, 3] = alpha * 255 * np.ones([256])
    L = len(rgb)
    N = 256 // L
    n = 0
    for i in range(L-1):
        color = rgb[i]
        alpha = opacity[i]
        lut[n:n+N, 0] = color[0]  * np.ones([N])
        lut[n:n+N, 1] = color[1]  * np.ones([N])
        lut[n:n+N, 2] = color[2]  * np.ones([N])
        lut[n:n+N, 3] = alpha * 255 * np.ones([N])
        n += N
    return lut

def bwr_lut_table(expansion=2):
    X = np.linspace(0, 1, 128)
    A = 255 * (X**expansion)
    B = A[::-1]
    C = 255 * np.ones(128)
    lut = np.ones((256,4)) * 255
    R = np.hstack((C,B))
    G = np.hstack((A,B))
    B = np.hstack((A,C))
    lut[:,0] = R
    lut[:,1] = G
    lut[:,2] = B
    return lut

def rgb_lut_table(expansion=2):
    X = np.linspace(0, 1, 64)
    increase = 255 * (X**(expansion**(-1)))
    decrease = increase[::-1]
    one = 255 * np.ones(64)
    zero = np.zeros(64)
    lut = np.ones((256, 4)) * 255
    B = np.hstack((one, decrease, zero, zero))
    G = np.hstack((increase, one, one, decrease))
    R = np.hstack((zero, zero, increase, one))
    lut[:,0] = R
    lut[:,1] = G
    lut[:,2] = B
    return lut

def ramp(values, slopes=1):
    counter = np.sign(slopes)
    slopes = abs(slopes)
    val = 0
    if counter < 0:
        val = 1
    vmax = np.max(values)
    vmin = np.min(values)
    start = vmin
    domain = vmax - vmin
    step = domain / slopes
    out = np.zeros(len(values))
    for i in range(slopes):
        mask = np.logical_and(values <= start+step, values >= start)
        out[mask] = counter*(values[mask] - start) / step + val
        counter *= -1
        start += step
        if val == 0:
            val = 1
        else:
            val = 0
    return out

def xyz_lut_table(points, variation=2):
    if variation == 0:
        a = b = c = 1
    else:
        a = variation + 1
        b = variation
        c = -variation
    x = 255 * ramp(points[:,0], a)
    y = 255 * ramp(points[:,1], -b)
    z = 255 * ramp(points[:,2], c)
    a = 255 * np.ones(len(points))
    lut = np.column_stack((x, y, z, a))
    return lut

#------------------------------------------------------------------------------
#                            Make the ModuleManager
#------------------------------------------------------------------------------

def make_module(color, opacity, lut_range, lut_expansion,
               reverse_lut, data_range, **kwargs):

    module = ModuleManager()
    if data_range is not None:
        module.scalar_lut_manager.use_default_range = False
        max_data = data_range[1]
        min_data = data_range[0]
        if max_data == min_data:
            max_data += 1
            min_data -= 1

        if type(lut_range) == str:

            if lut_range == '-:0:+':
                abs_data = max(abs(max_data),abs(min_data))
                lut_range = np.array([-abs_data, abs_data])

            elif lut_range == '0:+':
                lut_range = np.array([0, max_data])

            elif lut_range == '-:0':
                lut_range = np.array([min_data, 0])

            elif lut_range == '-:+':
                lut_range = np.array([min_data, max_data])

        module.scalar_lut_manager.data_range = lut_range

        if isinstance(color, np.ndarray):
            module.scalar_lut_manager.number_of_colors = len(color)
            module.scalar_lut_manager.lut.table = color

        elif color == 'bwr_e':
            lut_e = bwr_lut_table(lut_expansion)
            module.scalar_lut_manager.lut.table = lut_e

        elif color == 'rgb_e':
            lut_e = rgb_lut_table(lut_expansion)
            module.scalar_lut_manager.lut.table = lut_e

        elif color == 'xyz':
            points = kwargs.get('points')
            lut = xyz_lut_table(points)
            module.scalar_lut_manager.number_of_colors = len(points)
            module.scalar_lut_manager.lut.table = lut
            module.scalar_lut_manager.use_default_range = True

        elif type(color) == list or type(color) == tuple:
            lut_c = lut_table(color, opacity)
            module.scalar_lut_manager.lut.table = lut_c

        else:
            try:
                module.scalar_lut_manager.lut_mode = color
                module.scalar_lut_manager.reverse_lut = reverse_lut
            except:
                lut_c = lut_table(color, opacity)
                module.scalar_lut_manager.lut.table = lut_c
    else:
        lut_c = lut_table(color, opacity)
        module.scalar_lut_manager.lut.table = lut_c
        module.scalar_lut_manager.use_default_range = True

    return module

def make_lut_table(module, color, opacity, lut_expansion, reverse_lut, **kwargs):
    if color == 'bwr_e':
        lut_e = bwr_lut_table(lut_expansion)
        module.scalar_lut_manager.lut.table = lut_e

    elif color == 'rgb_e':
        lut_e = rgb_lut_table(lut_expansion)
        module.scalar_lut_manager.lut.table = lut_e

    elif type(color) == list or type(color) == tuple:
        lut_c = lut_table(color, opacity)
        module.scalar_lut_manager.lut.table = lut_c
        #print(dir(module.scalar_lut_manager))

    else:
        try:
            module.scalar_lut_manager.lut_mode = color
            module.scalar_lut_manager.reverse_lut = reverse_lut
            module.scalar_lut_manager._lut_mode_changed(color)
        except:
            lut_c = lut_table(color, opacity)
            module.scalar_lut_manager.lut.table = lut_c

def make_lut_range(lut_range, data_range):
    max_data = data_range[1]
    min_data = data_range[0]
    if max_data == min_data:
        max_data += 1
        min_data -= 1
    if type(lut_range) == str:
        if lut_range == '-:0:+':
            abs_data = max(abs(max_data),abs(min_data))
            lut_range = np.array([-abs_data, abs_data])
        elif lut_range == '0:+':
            lut_range = np.array([0, max_data])
        elif lut_range == '-:0':
            lut_range = np.array([min_data, 0])
        elif lut_range == '-:+':
            lut_range = np.array([min_data, max_data])
    else:
        lut_range = np.array(lut_range)
    return lut_range
#------------------------------------------------------------------------------
#                   Make the ModuleManager for Vectors
#------------------------------------------------------------------------------

def make_vector_module(color, opacity, lut_range, lut_expansion, reverse_lut):

    module = ModuleManager()

    if color == 'bwr_e':
        lut_e = bwr_lut_table(lut_expansion)
        module.vector_lut_manager.lut.table = lut_e

    elif color == 'rgb_e':
        lut_e = rgb_lut_table(lut_expansion)
        module.vector_lut_manager.lut.table = lut_e

    elif type(color) == list or type(color) == tuple:
        lut_c = lut_table(color,opacity)
        module.vector_lut_manager.lut.table = lut_c

    else:
        try:
            module.vector_lut_manager.lut_mode = color
            module.vector_lut_manager.reverse_lut = reverse_lut
        except:
            lut_c = lut_table(color, opacity)
            module.vector_lut_manager.lut.table = lut_c

    module.vector_lut_manager.use_default_range = True
    return module

#------------------------------------------------------------------------------
#https://matplotlib.org/1.2.1/mpl_examples/pylab_examples/show_colormaps.pdf
#------------------------------------------------------------------------------

def lut_mode_library():
    lut_mode = ['bwr_e','rgb_e', 'xyz',
    'Accent','Blues','Br_bG','Bu_gn','Bu_pu','CMRmap','Dark2','Gn_bu',
    'Greens','Greys','Or_rd','Oranges','PRGn','Paired','Pastel1','Pastel2',
    'Pi_yG','Pu_bu','Pu_bu_gn','Pu_or','Pu_rd','Purples','Rd_bu','Rd_gy','Rd_pu',
    'Rd_yl_bu','Rd_yl_gn','Reds','Set1','Set2','Set3','Spectral','Vega10','Vega20',
    'Vega20b','Vega20c','Wistia','Yl_gn','Yl_gn_bu','Yl_or_br','Yl_or_rd','afmhot',
    'autumn','binary','black-white','blue-red','bone','brg','bwr','cool',
    'coolwarm','copper','cubehelix','file','flag','gist_earth','gist_gray',
    'gist_heat','gist_ncar','gist_rainbow','gist_stern','gist_yarg','gnuplot',
    'gnuplot2','gray','hot','hsv','inferno','jet','magma','nipy_spectral',
    'ocean','pink','plasma','prism','rainbow','seismic','spectral','spring',
    'summer','terrain','viridis','winter',
    ]
    return lut_mode

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def check_arguments(**kwargs):
    point_data = kwargs.get('point_data')
    message = '*point_data* argument must be a list or a np.ndarray!'
    if point_data is not None:
        if type(point_data) is not np.ndarray and type(point_data) is not list:
            raise ValueError(message)
    vertex_data = kwargs.get('vertex_data')
    message = '*vertex_data* argument must be a list or a np.ndarray!'
    if vertex_data is not None:
        if type(vertex_data) is not np.ndarray and type(vertex_data) is not list:
            raise ValueError(message)
    face_data = kwargs.get('face_data')
    message = '*face_data* argument must be a list or a np.ndarray!'
    if face_data is not None:
        if type(face_data) is not np.ndarray and type(face_data) is not list:
            raise ValueError(message)
    edge_data = kwargs.get('edge_data')
    message = '*edge_data* argument must be a list or a np.ndarray!'
    if edge_data is not None:
        if type(edge_data) is not np.ndarray and type(edge_data) is not list:
            raise ValueError(message)
    indices = kwargs.get('vertex_indices')
    message = '*vertex_indices* argument must be a list or a np.ndarray!'
    if indices is not None:
        if type(indices) is not np.ndarray and type(indices) is not list:
            raise ValueError(message)
    color = kwargs.get('color')
    message = '*color* argument must be:\n'
    message += ' - an (r,g,b) tuple with [0,255] floats\n'
    message += ' - a color keyword string: '
    colors = color_library()
    message += str(colors.keys()) + '\n'
    message += ' - a palette keyword string: '
    palettes = palette_library()
    message += str(palettes.keys()) + '\n'
    message += ' - a list of (r,g,b) tuples and/or color keyword strings\n'
    message += ' - a color table keyword string: '
    tables = lut_mode_library()
    message += str(tables) + '\n'
    if color is not None:
        if isinstance(color, np.ndarray):
            pass
        else:
            if type(color) is not list:
                color = [color]
            for c in color:
                if type(c) is tuple and len(c) == 3:
                    for val in c:
                        if type(val) is float or type(val) is int:
                            if val >= 0 and val <= 255:
                                pass
                            else:
                                raise ValueError(message)
                        else:
                            raise ValueError(message)
                elif type(c) is str:
                    if c in colors:
                        pass
                    elif c in palettes:
                        pass
                    elif c in tables and len(color) == 1:
                        pass
                    else:
                        raise ValueError(message)
                else:
                    raise ValueError(message)
    opacity = kwargs.get('opacity')
    message = '*opacity* argument must be a [0,1] float!'
    if opacity is not None:
        if type(opacity) is float or type(opacity) is int:
            if opacity >= 0 and opacity <= 1:
                pass
            else:
                raise ValueError(message)
    glyph_type = kwargs.get('glyph_type')
    message = '*glyph_type* argument must be a keyword string: '
    types = ['cube', 'sphere', '3D-arrow', 'cylinder', 'line', 'cone', 'axes']
    message += str(types)
    if glyph_type is not None:
        if type(glyph_type) is str:
            if glyph_type in types:
                pass
            else:
                raise ValueError(message)
        else:
            raise ValueError(message)
    lut_range = kwargs.get('lut_range')
    message = '*lut_range* argument must be a list or tuple of length 2 '
    message += 'with float elements, '
    message += 'a float np.ndarray of shape (2,) or a keyword string: '
    ranges = ['-:+', '-:0:+', '-:0', '0:+']
    message += str(ranges)
    if lut_range is not None:
        if type(lut_range) is tuple or type(lut_range) is list:
            if len(lut_range) == 2:
                v1 = float(lut_range[0])
                v2 = float(lut_range[1])
                if type(v1) is float  and type(v2) is float :
                    pass
                else:
                    raise ValueError(message)
            else:
               raise ValueError(message)
        elif type(lut_range) is np.ndarray:
            print(lut_range.shape)
            if len( lut_range.shape) == 1 and lut_range.shape[0] == 2:
                pass
            else:
                raise ValueError(message)
        elif type(lut_range) is str:
            if lut_range in ranges:
                pass
            else:
                raise ValueError(message)
        else:
            raise ValueError(message)
    anchor_mode = kwargs.get('anchor_mode')
    message = '*anchor_mode* argument must be a keyword string: '
    types = ['vertex', 'face', 'edge']
    message += str(types)
    if anchor_mode is not None:
        if anchor_mode in types:
            pass
        else:
            raise ValueError(message)
