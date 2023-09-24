# -*- coding: utf-8 -*-
"""
Created on Thu May 26 10:20:13 2022

@author: WANGH0M
"""
__author__ = 'Hui'
#------------------------------------------------------------------------------
import numpy as np

from scipy import sparse
#------------------------------------------------------------------------------
from archgeolab.constraints.constraints_basic import column3D
# -------------------------------------------------------------------------

"""
from constraints.constraints_fairness import
con_fair_midpoint
con_laplacian_fairness
con_fairness_4th_different_polylines
"""

# -------------------------------------------------------------------------
#                           fairness
# -------------------------------------------------------------------------
def con_fair_midpoint(v,vl,vr,move,Vnum,N):
    "vl+vr-2v = 0"
    num = len(v)      
    arr = np.arange(3*num)
    one = np.ones(3*num)
    row = np.tile(arr,3)
    cc = column3D(v, move,Vnum)
    c1 = column3D(vl,move,Vnum)
    c2 = column3D(vr,move,Vnum)
    col = np.r_[cc,c1,c2]
    data = np.r_[2*one,-one,-one]
    K = sparse.coo_matrix((data,(row,col)), shape=(3*num, N))       
    return K

def con_laplacian_fairness(v,neib,move,Vnum,N):
    "v1+v2+..+vn = n*v"
    valence = neib.shape[1] # column num 
    num = len(v)   
    one = np.ones(3*num)
    row = np.tile(np.arange(3*num),valence+1)
    col = column3D(v, move,Vnum)
    data = valence*one
    for i in range(valence):
        ci = column3D(neib[:,i], move,Vnum)
        col = np.r_[col, ci]
        data = np.r_[data,-one]
    K = sparse.coo_matrix((data,(row,col)), shape=(3*num, N))       
    return K 

def con_fairness_4th_different_polylines(pylist,diag=False,**kwargs):
    "generate con_fairness_4th_different to any given polyline-list"
    "(v0-4*v1+6*v2-4*v3+v4)^2=0"
    if diag:
        w = kwargs.get('fairness_diag_4diff')
    else:
        w = kwargs.get('fairness_4diff')
    mesh = kwargs.get('mesh')
    X = kwargs.get('X')
    v0=v1=v2=v3=v4 = np.array([],dtype=int)
    for pl in pylist:
        if len(pl)>4:
            v0 = np.r_[v0,pl[:-4]]
            v1 = np.r_[v1,pl[1:-3]]
            v2 = np.r_[v2,pl[2:-2]]
            v3 = np.r_[v3,pl[3:-1]]
            v4 = np.r_[v4,pl[4:]]
    num = len(v0)
    arr = np.arange(num)
    row = np.tile(arr,15)
    c0 = column3D(v0,0,mesh.V)
    c1 = column3D(v1,0,mesh.V)
    c2 = column3D(v2,0,mesh.V)
    c3 = column3D(v3,0,mesh.V)
    c4 = column3D(v4,0,mesh.V)
    col = np.r_[c0,c1,c2,c3,c4]
    d = np.r_[X[c0]-4*X[c1]+6*X[c2]-4*X[c3]+X[c4]]
    data = 2*np.r_[d,-4*d,6*d,-4*d,d]
    r = np.linalg.norm(d.reshape(-1,3,order='F'),axis=1)**2
    H = sparse.coo_matrix((data,(row,col)), shape=(num, len(X)))
    return H*w,r*w
