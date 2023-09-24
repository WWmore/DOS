# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:55:43 2022

@author: WANGH0M
"""
__author__ = 'Hui'
#---------------------------------------------------------------------------
import numpy as np

from geometrylab.geometry import Polyline

#-----------------------------------------------------------------------------

"""from huilab.huimesh.curves import 
mesh_polylines,
make_polyline_from_endpoints,
make_multiple_polylines_from_endpoints,
get_isoline_between_2bdry,
get_diagonal_polyline_from_2points
"""
####-------------------------------------------------------------

    
def mesh_polylines(V, plys):
    "copy from meshpy.py/mesh_polylines; plys are seglists" #not work
    polylines = []
    for family in plys:
        poly_family = []
        for curve in family:
            poly_family.append(Polyline(V[curve,:]))
        polylines.append(poly_family)
    return polylines

def diagonal_cell_array(v1):
    """represent polyline data
    suppose num=6
    a0 = [2,2,2,2,2,2]
    a1 = [0,1,2,3,4,5]
    a2 = [6,7,8,9,10,11] 
    cells = [2,0,6, 2,1,7,..., 2,5,11]
    """
    num = len(v1)
    c = np.repeat(2,num)
    "i = np.arange(num), j = num + i"
    i,j = np.arange(2*num).reshape(2,-1)
    cells = np.vstack((c,i,j)).T
    cells = np.ravel(cells)
    return cells

def make_polyline_from_endpoints(Vl,Vr):
    "fork from quadring.py"
    VV = np.vstack((Vl,Vr))
    v = np.arange(len(Vl))
    data = diagonal_cell_array(v)
    poly = Polyline(VV)
    poly.cell_array=data
    #poly.refine(steps=3)  
    return poly

    
def get_diagonal_polyline_from_2points(mesh,vi,is_poly=True):
    if len(vi)==1:
        print('You should select at least two diagonal vertices!')
    else:
        H,V = mesh.halfedges, mesh.vertices
        vl = vr = np.array([],dtype=int)
        es1 = np.where([H[:,0]==vi[0]][0])[0]
        es2 = np.where([H[:,0]==vi[1]][0])[0]
        f1 = H[es1,1]
        f2 = H[es2,1]
        f = np.intersect1d(f1,f2)
        il = es1[np.where(f1==f)[0]]
        ir = es2[np.where(f2==f)[0]]

        vl = np.r_[vl,H[il,0]]
        num = len(np.where(H[:,0]==H[il,0])[0])
        while num==4 and H[H[H[H[il,4],2],4],1]!=-1 and H[H[il,4],1]!=-1:
            il = H[H[H[H[il,4],2],4],3]
            if H[il,0] in vl:
                break
            vl = np.r_[vl, H[il,0]]
            num = len(np.where(H[:,0]==H[il,0])[0])
        
        vr = np.r_[vr,H[ir,0]]
        num = len(np.where(H[:,0]==H[ir,0])[0])
        while num==4 and H[H[H[H[ir,4],2],4],1]!=-1 and H[H[ir,4],1]!=-1:
            ir =H[H[H[H[ir,4],2],4],3]
            if  H[ir,0] in vr:
                break
            vr = np.r_[vr, H[ir,0]] 
            num = len(np.where(H[:,0]==H[ir,0])[0])
        "Vl = self.vertices[vl[::-1]]; Vr = self.vertices[vr]"
        iv = np.r_[vl[::-1],vr]
        VV = V[iv]
        if is_poly:
            poly = make_polyline_from_endpoints(VV[:-1,:],VV[1:,:])
            return iv,VV,poly
        return iv,VV

    
def make_multiple_polylines_from_endpoints(VV,ns):
    def _multiple_segments_cell_array(ns):
        "ns is an array of nums of each segment"
        ci = np.array([],dtype=int)
        a = 0
        for n in ns:
            i = a + np.arange(n)
            ci = np.r_[ci,i]
            a += n + 1
        c = np.ones(len(ci))*2
        cj = ci + 1 
        cells = np.vstack((c,ci,cj)).T
        cells = np.ravel(cells)
        return cells          
    data = _multiple_segments_cell_array(ns)
    poly = Polyline(VV)
    poly.cell_array=data
    return poly   

def get_isoline_between_2bdry(mesh,vb):
    "starting from 1 boundary-vertex, get isoline untile to opposit bdry"
    "work for rectangular-patch + cylinder-annulus; vb in continious order"
    H = mesh.halfedges
    if True:
        "reverse the boundary list"
        vi,vj = vb[0],vb[1]
        ei = np.intersect1d(np.where(H[:,0]==vi)[0], np.where(H[H[:,4],0]==vj)[0])[0]
        if H[ei,1]!=-1:
            vb = vb[::-1]
    allv = []
    for i in range(len(vb)-1):
        vi,vj = vb[i],vb[i+1]
        ei = np.intersect1d(np.where(H[:,0]==vi)[0], np.where(H[H[:,4],0]==vj)[0])[0]
        plv = [vi]
        e = ei
        while H[H[e,4],1]!=-1:
            e = H[H[H[e,4],2],2]
            plv.append(H[e,0])
        allv.append(plv)  
    "last boundary"    
    plvj = [vj]    
    e = H[ei,4]
    while H[e,1]!=-1:
        e = H[H[H[e,2],2],4]
        plvj.append(H[e,0])
    allv.append(plvj)   
    return allv