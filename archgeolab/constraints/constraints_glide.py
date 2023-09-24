# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 11:29:34 2021

@author: WANGH0M
"""
__author__ = 'Hui'
#------------------------------------------------------------------------------
import numpy as np

from scipy import sparse
#------------------------------------------------------------------------------
from geometrylab import utilities

from archgeolab.constraints.constraints_basic import column3D
#------------------------------------------------------------------------------
"""
from constraints_glide import 
    con_alignment,con_alignments,con_glide_in_plane,
    con_sharp_corner,con_fix_vertices
    con_selected_vertices_glide_in_one_plane
"""
#------------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    #                             Glide / Alignment
    # -------------------------------------------------------------------------
def con_alignments(w,refPoly, glideInd, overlap=False,**kwargs):
    #w = kwargs.get('i_boundary_glide')
    N = kwargs.get('N')
    null = np.zeros([0])
    H = sparse.coo_matrix((null,(null,null)), shape=(0,N))
    r = np.array([])
    for i in range(len(glideInd)):
        v, crv = glideInd[i], refPoly[i]
        Hi,ri = con_alignment(w,crv,v,**kwargs)
        H = sparse.vstack((H,Hi))
        r = np.r_[r,ri]
    return H,r

def con_alignment(w, refPoly, glideInd, overlap=False,**kwargs):
    """glide on a given polylineCurve (may be not closed)
       refPoly: reference polylineCurve
       glideInd: indices of constraint vertices from the mesh
       constraints: from tangents e1 to get orthogonal e2,e3,
                    then(P-Q)*e2=0; (P-Q)*e3=0
                    such that P-Q align nealy on direction e1
            where P (varaibles), and Q,e2,e3 are computed_concrete_values
       linear equations: P*e2=Q*e2; P*e3=Q*e3
                ~       [P;P] * [e2;e3] = [Q;Q] * [e2;e3]
       another way: P-Q // e1 <==> (P-Q) x e1 = 0
    """
    #w = kwargs.get('boundary_glide')
    mesh = kwargs.get('mesh')
    X = kwargs.get('X')
    ind = glideInd
    Vr = refPoly.vertices
    Tr = refPoly.vertex_tangents()
    c_p = column3D(ind,0,mesh.V)
    P = X[c_p].reshape(-1,3,order='F')
    closeP = refPoly.closest_vertices(P)
    Q = Vr[closeP]
    e1 = Tr[closeP]
    e2 = utilities.orthogonal_vectors(e1)
    e3 = np.cross(e1,e2)
    r = np.r_[np.einsum('ij,ij->i',Q,e2),np.einsum('ij,ij->i',Q,e3)]

    num = len(ind)
    row = np.tile(np.arange(2*num),3)
    col = column3D(np.tile(ind,2),0,mesh.V)
    ee = np.vstack((e2,e3))
    data = ee.flatten('F')
    H = sparse.coo_matrix((data,(row,col)), shape=(2*num,len(X)))
    if overlap:
        "P=Q"
        row = np.arange(3*num)
        col = c_p
        data = np.ones(3*num)
        r0 = Q.flatten('F')
        H0 = sparse.coo_matrix((data,(row,col)), shape=(3*num,len(X)))
        H = sparse.vstack((H,H0*1))
        r = np.r_[r,r0*1]
    return H*w,r*w

def con_glide_in_plane(coo,**kwargs):
    "glide on xy plane: coo=0,1,2 ~ yz, xz, xy "
    w = kwargs.get('z0')
    mesh = kwargs.get('mesh')
    N = kwargs.get('N')
    numv = mesh.V
    row = np.arange(numv)
    crr = np.arange(numv)
    col = coo*numv+crr
    data = np.ones(numv)
    H = sparse.coo_matrix((data,(row,col)), shape=(numv, N))
    r = np.zeros(numv)
    return H*w, r*w

def con_selected_vertices_glide_in_one_plane(v,coo,value,**kwargs):
    "glide on x/y/z plane: coo=0,1,2 ~ yz, xz, xy, e.g. z=1 <==>coo=2,value=1"
    mesh = kwargs.get('mesh')
    N = kwargs.get('N')
    numv = len(v)
    row = np.arange(numv)
    col = coo*mesh.V+v
    data = np.ones(numv)
    H = sparse.coo_matrix((data,(row,col)), shape=(numv, N))
    r = np.ones(numv) * value
    return H, r

def con_sharp_corner(move=0,angle=90,**kwargs):
    "if given angle, need change (vl-v)*(vr-v) = 0"
    w = kwargs.get('sharp_corner')
    mesh = kwargs.get('mesh')
    X = kwargs.get('X')
    v = mesh.corner
    neib = []
    for i in v:
        neib.append(mesh.ringlist[i])
    neib = np.array(neib)
    vl,vr = neib[:,0],neib[:,1]
    c_v = column3D(v,move,mesh.V)
    c_vl = column3D(vl,move,mesh.V)
    c_vr = column3D(vr,move,mesh.V)
    col = np.r_[c_v,c_vl,c_vr]
    row = np.tile(np.arange(len(v)),9)
    data = np.r_[2*X[c_v]-X[c_vl]-X[c_vr],X[c_vr]-X[c_v],X[c_vl]-X[c_v]]
    r = np.einsum('ij,ij->i',(X[c_vl]-X[c_v]).reshape(-1,3, order='F'),(X[c_vr]-X[c_v]).reshape(-1,3, order='F'))  
    H = sparse.coo_matrix((data,(row,col)), shape=(len(v), len(X)))
    return H*w,r*w
    
def con_fix_vertices( index, Vf,**kwargs):
    "X[column(index)]==Vf"
    w = kwargs.get('fix_point')
    mesh = kwargs.get('mesh')
    N = kwargs.get('N')
    try:
        row = np.arange(3*len(index))
        num = 3*len(index)
    except:
        row = np.arange(3*1)
        num = 3
    col = column3D(index,0,mesh.V)
    data = np.ones(3*1)
    r = Vf
    H = sparse.coo_matrix((data,(row,col)), shape=(num,N))
    return H*w,r*w
    