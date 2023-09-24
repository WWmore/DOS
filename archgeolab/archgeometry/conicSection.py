# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:16:21 2022

@author: WANGH0M
"""
__author__ = 'Hui Wang'
# -----------------------------------------------------------------------------
import numpy as np

import scipy
#------------------------------------------------------------------------------
from geometrylab.geometry.meshprimitives import mesh_sphere
# -----------------------------------------------------------------------------
"""
vs_sphere_equation
sphere_equation
get_vs_interpolated_sphere
interpolate_sphere
get_sphere_packing
"""

#------------------------------------------------------------------------------
#                                 GENERATION
#------------------------------------------------------------------------------
def vs_sphere_equation(V,vneib):
    "V=vertices_coordinate, index: vneib=[v,neib]"
    "P = [x^2+y^2+z^2,x,y,z,1]"
    Pc0 = np.einsum('ij,ij->i',V[vneib],V[vneib])
    P = np.c_[Pc0, V[vneib], np.ones(len(vneib))]
    H = np.einsum('ij,jk',P.T,P)
    Q = np.array([[0,0,0,0,-2],[ 0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[-2,0,0,0,0]])
    vals, vecs = scipy.linalg.eig(H,Q)
    vals = np.abs(vals)
    smallest = list(vals).index(min(list(vals[vals>=0])))
    ####smallest = np.argmin(vals)
    vector = vecs[:,smallest]
    A, B, C, D, E = vector[0],vector[1],vector[2],vector[3],vector[4]
    delt = np.sqrt(np.abs(B*B+C*C+D*D-4*A*E))
    A,B,C,D,E = A/delt,B/delt,C/delt,D/delt,E/delt
    eps = np.finfo(float).eps
    A = A+eps
    cM = -1/(2*A) * np.array([B, C, D])
    r = np.abs((B*B+C*C+D*D-4*A*E)/(4*A*A))
    if A < 0:
        coo = [-A,-B,-C,-D,-E]
    else:
        coo = [A,B,C,D,E]
    # test: coo * V[vneib] * coo--->0
   # print 'test', vals[smallest], np.dot(np.dot(np.array(coo), H),np.array(coo).reshape(-1,1))
    return coo, cM, np.sqrt(r)

def sphere_equation(V,order,somelist):
    "somelist = ringlist / vringlist"
    V4 = V[order]
    coolist,clist,rlist = [],[],[]
    for v in order:
        ring = somelist[v]
        coo,c,r = vs_sphere_equation(V,ring)
        coolist.append(coo)
        clist.append(c)
        rlist.append(r)
    coo = np.array(coolist)
    sphere_coeff = (coo.T).flatten()
    M = np.array(clist).reshape(-1,3)
    nx = 2*coo[:,0]*V4[:,0]+coo[:,1]
    ny = 2*coo[:,0]*V4[:,1]+coo[:,2]
    nz = 2*coo[:,0]*V4[:,2]+coo[:,3]
    sphere_n = np.c_[nx,ny,nz].flatten()
    return M,np.array(rlist),sphere_coeff,sphere_n

def get_vs_interpolated_sphere(V,v,ringlist):
    C,r,_,_ = sphere_equation(V,v,ringlist)
    all_vi = np.array([],dtype=int)
    for i in v:
        all_vi = np.r_[all_vi,np.array(ringlist[i])]
    Vneib = V[all_vi]
    return C,r,Vneib
    
def interpolate_sphere(V0,V1,V2,V3,V4):
    """interpolate 5 vertices
    a(x^2+y^2+z^2)+(bx+cy+dz)+e=0 ; normalize: F^2 = b^2+c^2+d^2-4ae=1
    sphere center C:= (m1,m2,m3) = -(b, c, d) /a/2
    sphere radius:= F^2/a/2 = 1/(2a)
    shere normal N:// V-C // (Vx+B/A/2, Vy+C/A/2, Vz+D/A/2);
    unit_sphere_normal==-(2*A*Vx+B, 2*A*Vy+C, 2*A*Vz+D), (note direction: from vertex to center)
    since P(Vx,Vy,Vz) satisfy the sphere eq. and the normalizated eq., so that
        N ^2=1
    """
    num = len(V0)
    allV = np.vstack((V0,V1,V2,V3,V4))
    arr = num*np.arange(5)
    def _vs(partV):
        "P = [x^2+y^2+z^2,x,y,z,1]"
        Pc0 = np.einsum('ij,ij->i',partV,partV)
        P = np.c_[Pc0, partV, np.ones(5)]
        H = np.einsum('ij,jk',P.T,P)
        Q = np.array([[0,0,0,0,-2],[ 0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[-2,0,0,0,0]])
        vals, vecs = scipy.linalg.eig(H,Q)
        vals = np.abs(vals)
        smallest = list(vals).index(min(list(vals[vals>=0])))
        vector = vecs[:,smallest]
        A, B, C, D, E = vector[0],vector[1],vector[2],vector[3],vector[4]
        delt = np.sqrt(np.abs(B*B+C*C+D*D-4*A*E))
        A,B,C,D,E = A/delt,B/delt,C/delt,D/delt,E/delt
        cM = -1/(2*A) * np.array([B, C, D])
        r = np.abs((B*B+C*C+D*D-4*A*E)/(4*A*A))
        if A < 0:
            coo = [-A,-B,-C,-D,-E]
        else:
            coo = [A,B,C,D,E]
        return coo, cM, np.sqrt(r)

    coolist,clist,rlist = [],[],[]
    for i in range(num):
        partV = allV[i+arr]
        coo,c,r = _vs(partV)
        coolist.append(coo)
        clist.append(c)
        rlist.append(r)
    coo = np.array(coolist)
    sphere_coeff = (coo.T).flatten()
    M = np.array(clist).reshape(-1,3)
    nx = 2*coo[:,0]*V0[:,0]+coo[:,1] #==2*A*Vx+B
    ny = 2*coo[:,0]*V0[:,1]+coo[:,2] #==2*A*Vx+C
    nz = 2*coo[:,0]*V0[:,2]+coo[:,3] #==2*A*Vx+D
    sphere_n = -np.c_[nx,ny,nz]#.flatten() ##note the sign
    return M,np.array(rlist),sphere_coeff,sphere_n
       
def get_sphere_packing(C,r,Fa=20,Fv=20):
    num = C.shape[0]
    M0 = mesh_sphere(C[0],r[0],Fa,Fv)
    for i in range(num-1):
        Si = mesh_sphere(C[i+1], r[i+1],Fa,Fv)
        half = Si.halfedges
        V,E,F = Si.V, Si.E, Si.F
        M0.vertices = np.vstack((M0.vertices, Si.vertices))
        half[:,0] += (i+1)*V
        half[:,1] += (i+1)*F
        half[:,2] += (i+1)*2*E
        half[:,3] += (i+1)*2*E
        half[:,4] += (i+1)*2*E
        half[:,5] += (i+1)*2*E
        M0.halfedges = np.vstack((M0.halfedges, half))
    M0.topology_update()
    return M0

