# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:35:11 2021

@author: wangh0m
"""
__author__ = 'Hui'
#------------------------------------------------------------------------------
import numpy as np

from scipy import sparse
#------------------------------------------------------------------------------
"""
from constraints_basic import 
    column3D,con_edge,con_unit,con_bigger_than,\
    con_constl,con_equal_length,\
    con_planarity,con_planarity_constraints,con_unit_normal,con_orient
"""
# -------------------------------------------------------------------------
#                           general / basic
# -------------------------------------------------------------------------

def column3D(arr, num1, num2):
    """
    Parameters
    ----------
    array : array([1,4,7]).
    num1 : starting num.=100
    num2 : interval num.= 10

    Returns
    -------
    a : array(100+[1,4,7, 10,14,17, 20,24,27]).
    """
    a = num1 + np.r_[arr, num2+arr, 2*num2+arr]
    return a

def con_edge(X,c_v1,c_v3,c_ld1,c_ud1):
    "(v1-v3) = ld1*ud1"
    num = len(c_ld1)
    ld1 = X[c_ld1]
    ud1 = X[c_ud1]
    a3 = np.ones(3*num)
    row = np.tile(np.arange(3*num),4)
    col = np.r_[c_v1,c_v3,np.tile(c_ld1,3),c_ud1]
    data = np.r_[a3,-a3,-ud1,-np.tile(ld1,3)]
    r = -np.tile(ld1,3)*ud1
    H = sparse.coo_matrix((data,(row,col)), shape=(3*num, len(X)))
    return H,r

def con_unit_normal(X,c_e1,c_e2,c_e3,c_e4,c_n):
    "n^2=1; n*(e1-e3)=0; n*(e2-e4);"
    "Hui: better than (l*n=(e1-e3)x(e2-e4), but no orientation"
    H1,r1 = con_unit(X,c_n)
    H2,r2 = con_planarity(X,c_e1,c_e3,c_n)
    H3,r3 = con_planarity(X,c_e2,c_e4,c_n)
    H = sparse.vstack((H1,H2,H3))
    r = np.r_[r1,r2,r3]
    return H,r

def con_unit(X,c_ud1,w=100):
    "ud1**2=1"
    num = int(len(c_ud1)/3)
    arr = np.arange(num)
    row = np.tile(arr,3)
    col = c_ud1
    data = 2*X[col]
    r =  np.linalg.norm(X[col].reshape(-1,3,order='F'),axis=1)**2 + np.ones(num)
    H = sparse.coo_matrix((data,(row,col)), shape=(num, len(X)))
    return H*w,r*w

def con_constl(c_ld1,init_l1,N):
    "ld1 == const."
    num = len(c_ld1)
    row = np.arange(num,dtype=int)
    col = c_ld1
    data = np.ones(num,dtype=int)
    r = init_l1
    H = sparse.coo_matrix((data,(row,col)), shape=(num, N))
    return H,r

def con_bigger_than(X,minl,c_vi,c_vj,c_ai,num):
    "(vi-vj)^2-ai^2=minl"
    col = np.r_[c_vi,c_vj,c_ai]
    row = np.tile(np.arange(num),7)
    data = 2*np.r_[X[c_vi]-X[c_vj], -X[c_vi]+X[c_vj], -X[c_ai]]
    r = np.linalg.norm((X[c_vi]-X[c_vj]).reshape(-1,3,order='F'),axis=1)**2
    r = r - X[c_ai]**2 + np.ones(num)*minl
    H = sparse.coo_matrix((data,(row,col)), shape=(num, len(X)))
    return H,r

def con_planarity(X,c_v1,c_v2,c_n): 
    "n*(v1-v2)=0"
    num = int(len(c_n)/3)
    col = np.r_[c_n,c_v1,c_v2]
    row = np.tile(np.arange(num),9)
    data = np.r_[X[c_v1]-X[c_v2],X[c_n],-X[c_n]]
    r = np.einsum('ij,ij->i',X[c_n].reshape(-1,3, order='F'),(X[c_v1]-X[c_v2]).reshape(-1,3, order='F'))
    H = sparse.coo_matrix((data,(row,col)), shape=(num, len(X)))
    return H,r

def con_equal_length(X,c1,c2,c3,c4):
    "(v1-v3)^2=(v2-v4)^2"
    num = int(len(c1)/3)
    row = np.tile(np.arange(num),12)
    col = np.r_[c1,c2,c3,c4]
    data = 2*np.r_[X[c1]-X[c3],X[c4]-X[c2],X[c3]-X[c1],X[c2]-X[c4]]
    r = np.linalg.norm((X[c1]-X[c3]).reshape(-1,3, order='F'),axis=1)**2
    r = r-np.linalg.norm((X[c2]-X[c4]).reshape(-1,3, order='F'),axis=1)**2
    H = sparse.coo_matrix((data,(row,col)), shape=(num,len(X)))
    return H,r

def con_orient(X,Nv,c_vN,c_a,neg=False):
    "vN*Nv = a^2; if neg: vN*Nv = -a^2; variables: vN, a; Nv is given"
    if neg:
        sign = -1
    else:
        sign = 1    
    num = int(len(c_a))
    row = np.tile(np.arange(num),4)
    col = np.r_[c_vN,c_a]
    data = np.r_[Nv.flatten('F'),-sign*2*X[c_a]]
    r = -sign*X[c_a]**2
    H = sparse.coo_matrix((data,(row,col)), shape=(num,len(X))) 
    return H,r


    # -------------------------------------------------------------------------
    #                          Geometric Constraints (from Davide)
    # -------------------------------------------------------------------------

def con_normal_constraints(**kwargs):
    "represent unit normal: n^2=1"
    #w = kwargs.get('normal') * kwargs.get('geometric')
    mesh = kwargs.get('mesh')
    X = kwargs.get('X')
    V = mesh.V
    F = mesh.F
    f = 3*V + np.arange(F)
    i = np.arange(F)
    i = np.hstack((i, i, i)) #row ==np.tile(i,3) == np.r_[i,i,i]
    j = np.hstack((f, F+f, 2*F+f)) #col ==np.r_[f,F+f,2*F+f]
    data = 2 * np.hstack((X[f], X[F+f], X[2*F+f])) #* w
    H = sparse.coo_matrix((data,(i,j)), shape=(F,len(X)))
    r = ((X[f]**2 + X[F+f]**2 + X[2*F+f]**2) + 1) #* w
    return H,r


def con_planarity_constraints(is_unit_edge=False,**kwargs):
    "n*(vi-vj) = 0; Note: making sure normals is always next to V in X[V,N]"
    w = kwargs.get('planarity')
    mesh = kwargs.get('mesh')
    X = kwargs.get('X')
    V = mesh.V
    F = mesh.F
    f, v1, v2 = mesh.face_edge_vertices_iterators(order=True)
    if is_unit_edge:
        "f*(v1-v2)/length = 0, to avoid shrinkage of edges"
        num = len(v1)
        col_v1 = column3D(v1,0,V)
        col_v2 = column3D(v2,0,V)
        col_f = column3D(f,3*V,F)
        Ver = mesh.vertices
        edge_length = np.linalg.norm(Ver[v1]-Ver[v2],axis=1)
        row = np.tile(np.arange(num), 9)
        col = np.r_[col_f,col_v1,col_v2]
        l = np.tile(edge_length,3)
        data = np.r_[(X[col_v1]-X[col_v2])/l, X[col_f]/l, -X[col_f]/l]
        H = sparse.coo_matrix((data,(row, col)), shape=(num, len(X)))
        r = np.einsum('ij,ij->i',X[col_f].reshape(-1,3,order='F'),(X[col_v1]-X[col_v2]).reshape(-1,3,order='F'))
        r /= edge_length
        Hn,rn = con_unit(X,col_f,10*w)
        H = sparse.vstack((H*w,Hn))
        r = np.r_[r*w,rn]
    else:
        K = f.shape[0]
        f = 3*V + f
        r = ((X[v2] - X[v1]) * X[f] + (X[V+v2] - X[V+v1]) * X[F+f]
             + (X[2*V+v2] - X[2*V+v1]) * X[2*F+f] ) * w
        v1 = np.hstack((v1, V+v1, 2*V+v1))
        v2 = np.hstack((v2, V+v2, 2*V+v2))
        f = np.hstack((f, F+f, 2*F+f))
        i = np.arange(K)
        i = np.hstack((i, i, i, i, i, i, i, i, i))
        j = np.hstack((f, v2, v1))
        data = 2 * np.hstack((X[v2] - X[v1], X[f], -X[f])) * w
        H = sparse.coo_matrix((data,(i,j)), shape=(K, len(X)))
        Hn,rn = con_normal_constraints(**kwargs)
        H = sparse.vstack((H*w,Hn*w*10))
        r = np.r_[r*w,rn*w*10]
    return H,r
