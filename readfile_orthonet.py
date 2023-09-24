# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:16:21 2022

@author: WANGH0M
"""
#------------------------------------------------------------------------------
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import sys

path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(path)

#print(path)
#------------------------------------------------------------------------------



# -------------------------------------------------------------------------
#                              Run
# -------------------------------------------------------------------------
if __name__ == '__main__':

    a = path + r'\objs'

    pq = a + r'\obj_pq'
    anet = a + r'\obj_anet'
    snet = a + r'\obj_snet'
    equ = a +r'\obj_equilibrium'
    
    
    file = equ + r'\quad_dome.obj'

    #----------------------------------------

    '''Instantiate the sample component'''
    from opt_gui_orthonet import OrthoNet
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
