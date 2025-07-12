# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:16:21 2022

@author: WANGH0M
"""

# All rights reserved.
#
# This software is free for non-commercial, research and evaluation use 
# under the terms of the LICENSE.md file.
#
# For inquiries contact hui.wang.1@kaust.edu.sa

#------------------------------------------------------------------------------
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import sys

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(path)

#print(path) ##/Users/X/Github/DOS
#------------------------------------------------------------------------------



# -------------------------------------------------------------------------
#                              Run
# -------------------------------------------------------------------------
if __name__ == '__main__':
    ## if MacBook:
    #file = path + '/objs/obj_equilibrium' + '/quad_dome.obj'
    #path = '/users/wanghui/Desktop/geometrylab7/'

    a = path + r'/objs'

    pq = a + r'/obj_pq'
    anet = a + r'/obj_anet'
    # snet = a + r'/obj_snet'
    # equ = a +r'/obj_equilibrium'
    
    #file = pq + r'/conical1.obj' #heart.obj
    file = anet + r'/knet1.obj'
    #file = snet + r'/cmc1.obj'
    #file = equ + r'/quad_dome.obj'
    


    #----------------------------------------

    '''Instantiate the sample component'''
    from dos_3_opt import OrthoNet
    component = OrthoNet()

    '''Instantiate the main geolab application'''
    from archgeolab.archgeometry.gui_basic import GeolabGUI
    GUI = GeolabGUI()

    '''Add the component to geolab'''
    GUI.add_component(component)
    
    '''Open an obj file'''
    GUI.open_obj_file(file)
    
    '''Open another obj file'''
    #GUI.open_obj_file(reffile)
    
    '''Start geolab main loop'''
    GUI.start()

