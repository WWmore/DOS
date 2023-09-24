# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 19:46:48 2022

@author: WANGH0M
"""
__author__ = 'Hui'
#---------------------------------------------------------------------------
import numpy as np
#------------------------------------------------------------------------------
from geometrylab.geometry.meshpy import Mesh
#-----------------------------------------------------------------------------



def make_quad_mesh_pieces(P1,P2,P3,P4):
    "fork from quadring.py"
    vlist = np.vstack((P1,P2,P3,P4))
    num = len(P1)
    arr = np.arange(num)
    flist = np.c_[arr,arr+num,arr+2*num,arr+3*num].tolist()
    ck = Mesh()
    ck.make_mesh(vlist,flist)
    return ck