# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 09:40:14 2023

@author: WANGH0M
"""
#------------------------------------------------------------------------------
import numpy as np

from scipy import sparse

from geometrylab import utilities
#------------------------------------------------------------------------------

"""
below constraints are copied from geometrylab/optimization/tuidedprojection.py
normal_constraints
planarity_constraints

edge_length_constraints
equilibrium_constraints,
compression_constraints,
area_constraints,
vector_area_constraints,
circularity_constraints

boundary_densities_constraints
fixed_boundary_normals_constraints
"""

# -------------------------------------------------------------------------
#                          Principal Meshes
# -------------------------------------------------------------------------

def normal_constraints(**kwargs): ## replaced by Hui's constraint
    "note: based on self.planarity=True: face_normal^2==1 "
    w = kwargs.get('normal') * kwargs.get('geometric')
    mesh = kwargs.get('mesh')
    V = mesh.V
    F = mesh.F
    N = kwargs.get('N')
    f = 3*V + np.arange(F)
    i = np.arange(F)
    i = np.hstack((i, i, i))
    j = np.hstack((f, F+f, 2*F+f))
    X = kwargs.get('X')
    data = 2 * np.hstack((X[f], X[F+f], X[2*F+f])) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(F, N))
    r = ((X[f]**2 + X[F+f]**2 + X[2*F+f]**2) + 1) * w
    #self.add_iterative_constraint(H, r, 'face_normal_length')
    print('n:', np.sum(np.square((H*X)-r)))
    return H,r

def planarity_constraints(**kwargs): ## replaced by Hui's constraint
    "note:  based on self.planarity=True: face_normal _|_ (vi-vj)"
    w = kwargs.get('planarity')
    mesh = kwargs.get('mesh')
    V = mesh.V
    F = mesh.F
    N = kwargs.get('N')
    X = kwargs.get('X')
    f, v1, v2 = mesh.face_edge_vertices_iterators(order=True)
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
    H = sparse.coo_matrix((data,(i,j)), shape=(K, N))
    #self.add_iterative_constraint(H, r, 'face_normal (planarity)')
    print('pq:', np.sum(np.square((H*X)-r)))
    return H,r

def circularity_constraints(**kwargs):
    "note:  based on self.circularity=True (N4)"
    w = kwargs.get('circularity')
    mesh = kwargs.get('mesh')
    V = mesh.V
    F = mesh.F
    N = kwargs.get('N')
    N3 = kwargs.get('N3')
    X = kwargs.get('X')
    f, v1, v2 = mesh.face_edge_vertices_iterators(order=True)
    K = f.shape[0]
    cx = N3 + np.array(f)
    cy = N3 + np.array(f) + F
    cz = N3 + np.array(f) + 2*F
    v1x = np.array(v1)
    v1y = np.array(V + v1)
    v1z = np.array(2*V + v1)
    v2x = np.array(v2)
    v2y = np.array(V + v2)
    v2z = np.array(2*V + v2)
    jx = np.hstack((v1x, v2x, cx))
    jy = np.hstack((v1y, v2y, cy))
    jz = np.hstack((v1z, v2z, cz))
    datax = np.hstack((X[v1x] - X[cx], -X[v2x] + X[cx], -X[v1x] + X[v2x]))
    rx = 0.5*X[v1x]**2 - 0.5*X[v2x]**2 -X[v1x]*X[cx] + X[v2x]*X[cx]
    datay = np.hstack((X[v1y] - X[cy], -X[v2y] + X[cy], -X[v1y] + X[v2y]))
    ry = 0.5*X[v1y]**2 - 0.5*X[v2y]**2 -X[v1y]*X[cy] + X[v2y]*X[cy]
    dataz = np.hstack((X[v1z] - X[cz], -X[v2z] + X[cz], -X[v1z] + X[v2z]))
    rz = 0.5*X[v1z]**2 - 0.5*X[v2z]**2 -X[v1z]*X[cz] + X[v2z]*X[cz]
    i = np.arange(K)
    i = np.hstack((i,i,i,i,i,i,i,i,i))
    j = np.hstack((jx, jy, jz))
    data = np.hstack((datax, datay, dataz)) * w
    r = (rx + ry + rz) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(K, N))
    #self.add_iterative_constraint(H, r,'circularity')
    return H,r

# -------------------------------------------------------------------------
#                             Equilibrium
# default these three run
# equilibrium_constraints; edge_length_constraints; boundary_densities_constraints
# -------------------------------------------------------------------------

def equilibrium_constraints(**kwargs):
    """based on self.equilibrium=True (N2); X[N1: N1+3E]
        edge_length:= X[N1: N1+E]; 
        density wij:= X[N1+E: N1+2E]
        sqrt_density sqrt(wij):= X[N1+2E: N1+3E]
    constraints: wij * (vi-vj) = [0; 0; 0.5*edge_length*Fe-P[:,2]]
    where edge_length is defined(constrained) in edge_length_constraints()
          density is defined(constrained) here equilibrium_constraints()
          sqrt_density is defined(constrained) in compression_constraints() (may not be used)
    """
    w = kwargs.get('equilibrium')
    mesh = kwargs.get('mesh')
    P = mesh.applied_forces
    Fe, Fa = mesh.self_weigth_loads ##note: known/computed values
    #print(P,Fe,Fa)
    norm = max(np.max(np.linalg.norm(P, axis=1)), np.max(Fe), np.max(Fa))
    w = 0.1 * w / (norm + 1e-6)
    V = mesh.V
    E = mesh.E
    F = mesh.F
    N = kwargs.get('N')
    N1 = kwargs.get('N1')
    X = kwargs.get('X')
    constrained = mesh.constrained_vertices
    inner = np.invert(np.in1d(np.arange(V), constrained))
    v0, vj = mesh.vertex_ring_vertices_iterators(sort=True)
    __, ej = mesh.vertex_ring_edges_iterators(sort=True)
    v_mask = np.invert(np.in1d(v0, constrained))
    v0 = v0[v_mask]
    vj = vj[v_mask]
    ej = ej[v_mask]
    ix = np.hstack((v0, v0, v0))
    jx = np.hstack((v0, vj, N1+E+ej))
    datax = np.hstack((X[N1+E+ej], -X[N1+E+ej], X[v0]-X[vj]))
    ri = X[N1+E+ej]*(X[v0]-X[vj])
    ri = utilities.sum_repeated(ri,v0)
    rx = -P[:,0]
    rx[inner] += ri
    iy = ix + V
    jy = np.hstack((V+v0, V+vj, N1+E+ej))
    datay = np.hstack((X[N1+E+ej], -X[N1+E+ej], X[V+v0]-X[V+vj]))
    ri = X[N1+E+ej]*(X[V+v0]-X[V+vj])
    ri = utilities.sum_repeated(ri,v0)
    ry = -P[:,1]
    ry[inner] += ri
    iz = np.hstack((iy+V, 2*V+v0))
    jz = np.hstack((2*V+v0, 2*V+vj, N1+E+ej, N1+ej)) ##desity,edge_length
    dataz = np.hstack((X[N1+E+ej], -X[N1+E+ej],
                       X[2*V+v0]-X[2*V+vj], -0.5*Fe[ej]))
    ri = X[N1+E+ej]*(X[2*V+v0]-X[2*V+vj])
    ri = utilities.sum_repeated(ri,v0)
    rz = -P[:,2]
    rz[inner] += ri
    i = np.hstack((ix, iy, iz))
    j = np.hstack((jx, jy, jz))
    data = np.hstack((datax, datay, dataz)) * w
    if kwargs.get('area') != 0:
        "note:  based on self.area=True (3)"
        N2 = kwargs.get('N2')
        vf, fj = mesh.vertex_ring_faces_iterators(sort=True)
        f_mask = np.invert(np.in1d(vf, constrained))
        vf = vf[f_mask]
        fj = fj[f_mask]
        Lf = mesh.face_lengths()
        iz = 2*V + vf
        jz = N2 + 3*F + fj
        dataz =  - Fa[fj] / Lf[fj] * w
        i = np.hstack((i,iz))
        j = np.hstack((j,jz))
        data = np.hstack((data, dataz))
    r = np.hstack((rx, ry, rz)) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(3*V, N))
    #self.add_iterative_constraint(H, r,'equilibrium')
    #print('eq:', np.sum(np.square((H*X)-r)))
    return H,r

def edge_length_constraints(**kwargs):
    "if self.equilibrium=True (N2); then it's used!"
    "variables edge_length E; constraints: E^2:= (Vi-Vj)^2;"
    w = kwargs.get('edge_length') * kwargs.get('geometric')
    mesh = kwargs.get('mesh')
    V = mesh.V
    E = mesh.E
    N = kwargs.get('N')
    N1 = kwargs.get('N1')
    X = kwargs.get('X')
    i = np.arange(E)
    e = N1 + i
    v1, v2 = mesh.edge_vertices()
    r = ( X[v1]**2 + X[v2]**2 - 2*X[v1]*X[v2]
        + X[V+v1]**2 + X[V+v2]**2 - 2*X[V+v1]*X[V+v2]
        + X[2*V+v1]**2 + X[2*V+v2]**2 - 2*X[2*V+v1]*X[2*V+v2]
        - X[e]**2 ) * w
    v1 = np.hstack((v1,V+v1,2*V+v1))
    v2 = np.hstack((v2,V+v2,2*V+v2))
    i = np.hstack((i,i,i,i,i,i,i))
    j = np.hstack((v1,v2,e))
    data = 2 * np.hstack((X[v1] - X[v2], X[v2] - X[v1], -X[e])) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(E, N))
    #self.add_iterative_constraint(H, r, 'edge_length')
    #print('el:', np.sum(np.square((H*X)-r)))
    return H,r

def boundary_densities_constraints(**kwargs):
    """if self.equilibrium=True (N2); then it's used!
       density wij:= X[N1+E: N1+2E];
       where boundary density=0; constraint: bdry_desity^2=0
    """
    w = kwargs.get('equilibrium')
    mesh = kwargs.get('mesh')
    E = mesh.E
    N = kwargs.get('N')
    N1 = kwargs.get('N1')
    X = kwargs.get('X')
    O = N1 + E
    constrained = mesh.constrained_vertices
    v1, v2 = mesh.edge_vertices()
    mask1 = np.in1d(v1, constrained)
    mask2 = np.in1d(v2, constrained)
    mask = np.logical_and(mask1, mask2)
    #mask = mesh.are_boundary_edges()
    j = O + np.arange(E)[mask]
    W = len(j)
    i = np.arange(len(j))
    data = 2*X[j] * w
    r = X[j]**2 * w
    H = sparse.coo_matrix((data,(i,j)), shape=(W,N))
    #self.add_constant_constraint(H, r, 'boundary_equilibrium')
    #print('bdry:', np.sum(np.square((H*X)-r)))
    return H,r

def compression_constraints(w, **kwargs):
    """if self.equilibrium=True (N2); it may use or not
        density wij:= X[N1+E: N1+2E]
        sqrt_density sqrt(wij):= X[N1+2E: N1+3E]
        constraint: density  = sqrt_density ^2
    """
    mesh = kwargs.get('mesh')
    sign = np.sign(w)
    E = mesh.E
    N = kwargs.get('N')
    N1 = kwargs.get('N1')
    X = kwargs.get('X')
    e = np.arange(E)
    i = np.hstack((e,e))
    j = np.hstack((N1+2*E+e, N1+E+e))
    data = np.hstack((2.0 * X[N1+2*E+e], -sign*np.ones(E))) * abs(w)
    H = sparse.coo_matrix((data,(i,j)), shape=(E,N))
    r = X[N1+2*E+e]**2 * abs(w)
    #self.add_iterative_constraint(H, r, 'compression')
    print('sqrt:', np.sum(np.square((H*X)-r)))
    return H, r


##below area from X[N3: N3+4F] are not used as mentioned in Davide's AAG-paper
def area_constraints(**kwargs):
    "note:  based on self.area=True (N3)"
    "Ax,Ay,Az := v1 x v2"
    w = kwargs.get('area') * kwargs.get('geometric')
    mesh = kwargs.get('mesh')
    V = mesh.V
    F = mesh.F
    N = kwargs.get('N')
    N2 = kwargs.get('N2')
    X = kwargs.get('X')
    fi, v1, v2 = mesh.face_edge_vertices_iterators(order=True)
    v1a = np.hstack((V+v1,2*V+v1,v1))
    v1b = np.hstack((2*V+v1,v1,V+v1))
    v2a = np.hstack((V+v2,2*V+v2,v2))
    v2b = np.hstack((2*V+v2,v2,V+v2))
    r = 0.5 * (X[v1a] * X[v2b] - X[v1b] * X[v2a]) * w
    k = np.hstack((fi,F+fi,2*F+fi))
    r = utilities.sum_repeated(r,k)
    f = N2 + np.arange(F)
    i = np.hstack((k,k,k,k,np.arange(3*F)))
    j = np.hstack((v1a, v2b, v1b, v2a, f, F+f, 2*F+f))
    data1 = 0.5 * np.hstack((X[v2b], X[v1a], -X[v2a], -X[v1b]))
    data2 = - np.ones(3*F)
    data = np.hstack((data1, data2)) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(3*F, N))
    #self.add_iterative_constraint(H, r, 'face_vector_area')
    print('a1:', np.sum(np.square((H*X)-r)))
    return H,r

def vector_area_constraints(**kwargs):
    "note:  based on self.equilibrium=True (N2)"
    "area^2:= Ax^2 + Ay^2 + Az^2"
    w = kwargs.get('area') * kwargs.get('geometric')
    mesh = kwargs.get('mesh')
    F = mesh.F
    N = kwargs.get('N')
    N2 = kwargs.get('N2')
    f = N2 + np.arange(F)
    i = np.arange(F)
    i = np.hstack((i,i,i,i))
    j = np.hstack((f,F+f,2*F+f,3*F+f))
    X = kwargs.get('X')
    data = 2 * np.hstack((X[f], X[F+f], X[2*F+f], - X[3*F+f])) * w
    H = sparse.coo_matrix((data,(i,j)), shape=(F, N))
    r = (X[f]**2 + X[F+f]**2 + X[2*F+f]**2 - X[3*F+f]**2) * w
    #self.add_iterative_constraint(H, r, 'face_area')
    print('a2:', np.sum(np.square((H*X)-r)))
    return H,r

# -------------------------------------------------------------------------
#                                 Boundary
# -------------------------------------------------------------------------

def fixed_boundary_normals_constraints(**kwargs): #not necessary, may not be used
    "based on self.planarity=True (N1)"
    "face_normal [Nx,Ny,Nz] from X[3V:3V+3F]; constraint: [Nx,Ny,Nz] are fixed"
    w = kwargs.get('fixed_boundary_normals')
    mesh = kwargs.get('mesh')
    N = kwargs.get('N')
    F = mesh.F
    V = mesh.V
    f = mesh.boundary_faces()
    X = kwargs.get('X')
    W = 3*len(f)
    nx = 3*V + f
    ny = 3*V + F + f
    nz = 3*V + 2*F + f
    data = np.ones(W) * w
    i = np.arange(W)
    j = np.hstack((nx, ny, nz))
    H = sparse.coo_matrix((data,(i,j)), shape=(W,N))
    r = np.hstack((X[nx],X[ny],X[nz])) * w
    #self.add_constant_constraint(H, r, 'fixed_boundary_normals')
    print('fix:', np.sum(np.square((H*X)-r)))
    return H, r
