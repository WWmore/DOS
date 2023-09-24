#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

import numpy as np

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.geometry.meshpy import Mesh

# -----------------------------------------------------------------------------

'''_'''

__author__ = 'Davide Pellis'

#------------------------------------------------------------------------------
#                                 GENERATION
#------------------------------------------------------------------------------

def mesh_plane(Fx=10, Fy=10, x_min=-5, x_max=5, y_min=-5, y_max=5):
    Fx +=1; Fy +=1
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M,(Fx,Fy))
    Q = np.zeros((Fx-1, Fy-1,4) ,dtype=np.int)
    Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,2] = M[:M.shape[0]-1,1:]
    Q[:,:,1] = M[1:,1:]
    Q[:,:,0] = M[1:,:M.shape[1]-1]
    Q = Q.reshape(((Fx-1)*(Fy-1),4))
    x = np.linspace(x_min,x_max,Fx)
    y = np.linspace(y_min,y_max,Fy)
    Px, Py = np.meshgrid(x,y)
    Px = np.reshape(Px,(Fx*Fy), order='F')
    Py = np.reshape(Py,(Fx*Fy), order='F')
    P = np.vstack((Px,Py,np.zeros(Fx*Fy))).T
    grid = Mesh()
    grid.make_mesh(P,Q)
    return grid

def mesh_cylinder(C=[0,0,0], r=1, h=3, Fa=10, Fv=5):
    C = np.array(C)
    Fx = Fa; Fy = Fv+1
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M,(Fy,Fx))
    M = np.insert(M,Fx,M[:,0],axis=1)
    Q = np.zeros((Fy-1, Fx,4) ,dtype=np.int)
    Q[:,:,0] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,1] = M[:M.shape[0]-1,1:]
    Q[:,:,2] = M[1:,1:]
    Q[:,:,3] = M[1:,:M.shape[1]-1]
    Q = Q.reshape(((Fx)*(Fy-1),4),order='C')
    phi = np.linspace(0,2*np.pi,Fx+1)
    phi = np.delete(phi,-1)
    z = np.linspace(C[2],C[2]+h,Fy)
    phi, z = np.meshgrid(phi,z)
    Px = np.cos(phi) + C[0]
    Py = np.sin(phi) + C[1]
    Px = np.reshape(Px,(Fx*(Fy)), order='F')
    Py = np.reshape(Py,(Fx*(Fy)), order='C')
    Pz = np.reshape(z,(Fx*(Fy)), order='C')
    P = np.vstack((Px,Py,Pz)).T
    grid = Mesh()
    grid.make_mesh(P,Q)
    return grid

def mesh_sphere(C=[0,0,0], r=1, Fa=20, Fv=10):
    C = np.array(C)
    Fx = Fa; Fy = Fv-1
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M, (Fy, Fx)) + 1
    M = np.insert(M, Fx, M[:,0], axis=1)
    Q = np.zeros((Fy-1, Fx, 4) ,dtype=np.int)
    Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,2] = M[:M.shape[0]-1,1:]
    Q[:,:,1] = M[1:,1:]
    Q[:,:,0] = M[1:,:M.shape[1]-1]
    Q = Q.reshape(((Fx)*(Fy-1),4),order='C').tolist()
    T = np.zeros((2, Fx, 3) ,dtype=np.int)
    T[0,:,0] = M[0,:M.shape[1]-1]
    T[0,:,1] = M[0,1:]
    T[0,:,2] = 0
    T[1,:,0] = M[-1,1:]
    T[1,:,1] = M[-1,:M.shape[1]-1]
    T[1,:,2] = Fx*Fy + 1
    T = T.reshape((Fx*2,3),order='C').tolist()
    F = Q + T
    phi = np.linspace(0,2*np.pi,Fx+1)[:-1]
    theta = np.linspace(0,np.pi,Fy+2)
    theta = theta[1:theta.shape[0]-1]
    phi, theta = np.meshgrid(phi,theta)
    Px = r*np.cos(phi)*np.sin(theta) + C[0]
    Py = r*np.sin(phi)*np.sin(theta) + C[1]
    Pz = r*np.cos(theta) + C[2]
    Px = np.reshape(Px,(Fx*(Fy)), order='C')
    Py = np.reshape(Py,(Fx*(Fy)), order='C')
    Pz = np.reshape(Pz,(Fx*(Fy)), order='C')
    P = np.vstack((Px,Py,Pz)).T
    P = np.insert(P,0,[C[0],C[1],C[2]+r], axis=0)
    P = np.insert(P,P.shape[0],[C[0],C[1],C[2]-r], axis=0)
    grid = Mesh()
    grid.make_mesh(P,F)
    return grid

def mesh_torus(C=[0,0,0], r1=3, r2=1, F1=25, F2=15):
    C = np.array(C)
    Fx = F2; Fy = F1
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M,(Fy,Fx))
    M = np.insert(M,Fx,M[:,0],axis=1)
    M = np.insert(M,Fy,M[0,:],axis=0)
    Q = np.zeros((Fy, Fx,4) ,dtype=np.int)
    Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,2] = M[:M.shape[0]-1,1:]
    Q[:,:,1] = M[1:,1:]
    Q[:,:,0] = M[1:,:M.shape[1]-1]
    Q = Q.reshape(((Fx)*(Fy),4),order='C')
    u = np.linspace(0, 2*np.pi, Fx+1)[:-1]
    v = np.linspace(0, 2*np.pi, Fy+1)[:-1]
    U, V = np.meshgrid(u, v)
    P = np.array([r2*np.cos(U)*np.cos(V) + r1*np.cos(V) + C[0],
                  r2*np.cos(U)*np.sin(V) + r1*np.sin(V) + C[1],
                  r2*np.sin(U) + C[2]])
    Px = np.reshape(P[0],(Fx*Fy), order='C')
    Py = np.reshape(P[1],(Fx*Fy), order='C')
    Pz = np.reshape(P[2],(Fx*Fy), order='C')
    P = np.vstack((Px,Py,Pz)).T
    grid = Mesh()
    grid.make_mesh(P,Q)
    return grid

def mesh_arrow(Fa=15, tip_length=3, shaft_length=8, tip_radius=1, shaft_radius=0.4):
    Fx = Fa; Fy = 7
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M,(Fy,Fx))
    M = np.insert(M,Fx,M[:,0],axis=1)
    Q = np.zeros((Fy-1, Fx,4) ,dtype=np.int)
    Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,2] = M[:M.shape[0]-1,1:]
    Q[:,:,1] = M[1:,1:]
    Q[:,:,0] = M[1:,:M.shape[1]-1]
    T = np.zeros((1, Fx, 3) ,dtype=np.int)
    T[0,:,0] = M[-1,1:]
    T[0,:,1] = M[-1,:M.shape[1]-1]
    T[0,:,2] = Fx*Fy
    Q = Q.reshape(((Fx)*(Fy-1),4),order='C').tolist()
    T = T.reshape((Fx,3),order='C').tolist()
    F = Q + T
    phi = np.linspace(0, 2*np.pi, Fa+1)[:-1]
    P0 = np.zeros((Fa,3))
    P1 = np.array([tip_radius*np.sin(phi), tip_radius*np.cos(phi),
                   -tip_length*np.ones(Fa)]).T
    P2 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                   -tip_length*np.ones(Fa)]).T
    P5 = np.array([0.3*tip_radius*np.sin(phi),0.3*tip_radius*np.cos(phi),
                   -0.3*tip_length*np.ones(Fa)]).T
    P6 = np.array([0.6*tip_radius*np.sin(phi),0.6*tip_radius*np.cos(phi),
                   -0.6*tip_length*np.ones(Fa)]).T
    P7 = np.array([0.1*tip_radius*np.sin(phi),0.1*tip_radius*np.cos(phi),
                   -0.1*tip_length*np.ones(Fa)]).T
    P3 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                   -shaft_length*np.ones(Fa)]).T
    P4 = np.array([[0,0,-shaft_length]])
    P = np.vstack((P0,P7,P5,P6,P1,P2,P3,P4))
    grid = Mesh()
    grid.make_mesh(P,F)
    return grid

def mesh_arrows(vectors, anchors, Fa=3, anchor_mode='head',
               tip_length=3, shaft_length=8, tip_radius=1, shaft_radius=0.4):
    Fx = Fa; Fy = 7
    if anchor_mode == 'tail':
        anchors = anchors + vectors
    L = np.linalg.norm(vectors,axis=1)
    N = vectors.shape[0]
    M = np.arange(Fx*Fy, dtype=np.int)
    M = np.reshape(M,(Fy,Fx))
    M = np.insert(M,Fx,M[:,0],axis=1)
    Q = np.zeros((Fy-1, Fx,4) ,dtype=np.int)
    Q[:,:,3] = M[:M.shape[0]-1,:M.shape[1]-1]
    Q[:,:,2] = M[:M.shape[0]-1,1:]
    Q[:,:,1] = M[1:,1:]
    Q[:,:,0] = M[1:,:M.shape[1]-1]
    T = np.zeros((1, Fx, 3) ,dtype=np.int)
    T[0,:,0] = M[-1,1:]
    T[0,:,1] = M[-1,:M.shape[1]-1]
    T[0,:,2] = Fx*Fy
    O = np.arange(0,N*Fa*Fy,Fa*Fy+1)
    O = np.reshape(O,(N,1,1))
    Q = Q.reshape(((Fx)*(Fy-1),4),order='C')
    T = T.reshape((Fx,3),order='C')
    Qd = np.zeros((N,Q.shape[0],Q.shape[1]),dtype=np.int)
    Qd[:] = Q
    Qd += O
    Td = np.zeros((N,T.shape[0],T.shape[1]),dtype=np.int)
    Td[:] = T
    Td += O
    Qd = np.insert(Qd,0,4,axis=2)
    Td = np.insert(Td,0,3,axis=2)
    Qd = np.reshape(Qd,(Qd.shape[0]*Qd.shape[1]*Qd.shape[2]))
    Td = np.reshape(Td,(Td.shape[0]*Td.shape[1]*Td.shape[2]))
    cells = np.hstack((Qd,Td))
    phi = np.linspace(0, 2*np.pi, Fa+1)[:-1]
    P0 = np.zeros((Fa,3))
    P1 = np.array([tip_radius*np.sin(phi), tip_radius*np.cos(phi),
                   -tip_length*np.ones(Fa)]).T
    P2 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                   -tip_length*np.ones(Fa)]).T
    P5 = np.array([0.3*tip_radius*np.sin(phi),0.3*tip_radius*np.cos(phi),
                   -0.3*tip_length*np.ones(Fa)]).T
    P6 = np.array([0.6*tip_radius*np.sin(phi),0.6*tip_radius*np.cos(phi),
                   -0.6*tip_length*np.ones(Fa)]).T
    P7 = np.array([0.1*tip_radius*np.sin(phi),0.1*tip_radius*np.cos(phi),
                   -0.1*tip_length*np.ones(Fa)]).T
    P3 = np.array([shaft_radius*np.sin(phi), shaft_radius*np.cos(phi),
                   -np.ones(Fa)]).T
    P4 = np.array([[0,0,-1]])
    P = np.vstack((P0,P7,P5,P6,P1,P2,P3,P4))
    Pt = np.zeros((N,P.shape[0],3))
    Pt[:] = P
    Off = np.tile(anchors,P.shape[0])
    Off = np.reshape(Off,(P.shape[0]*N,3))
    Ld = np.repeat(L,Fa+1)
    Ld = np.reshape(Ld,(N,Fa+1))
    Pt[:,6*Fa:,2] *= Ld
    Z = vectors/np.linalg.norm(vectors,axis=1,keepdims=True)
    X = utilities.orthogonal_vectors(Z)
    Y = np.cross(Z,X,axis=1)
    J = np.zeros((N,3,3))
    J[:,0,:] = X
    J[:,1,:] = Y
    J[:,2,:] = Z
    Pt = np.einsum('ijk,ikl->ijl',Pt,J)
    Pt = np.reshape(Pt,(P.shape[0]*N,3))
    Pt = Pt + Off
    return Pt, cells