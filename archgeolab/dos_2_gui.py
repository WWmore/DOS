# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:16:21 2022

@author: WANGH0M
"""
__author__ = 'Hui'
#------------------------------------------------------------------------------
import numpy as np

# -----------------------------------------------------------------------------
from geometrylab.optimization.guidedprojectionbase import GuidedProjectionBase

from archgeolab.constraints.constraints_basic import con_planarity_constraints

from archgeolab.constraints.constraints_fairness import con_fairness_4th_different_polylines

from archgeolab.constraints.constraints_net import con_unit_edge,\
    con_orthogonal_midline,con_anet,con_anet_diagnet,con_snet,\
    con_snet_diagnet,con_multinets_orthogonal

from archgeolab.constraints.constraints_glide import con_glide_in_plane,\
    con_alignment,con_alignments,con_selected_vertices_glide_in_one_plane,\
    con_fix_vertices,con_sharp_corner
            
from archgeolab.constraints.constraints_equilibrium import edge_length_constraints,\
    equilibrium_constraints,compression_constraints,area_constraints,\
    vector_area_constraints,\
    boundary_densities_constraints,fixed_boundary_normals_constraints
                    
from archgeolab.archgeometry.conicSection import interpolate_sphere

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

class GP_OrthoNet(GuidedProjectionBase):
    _N1 = 0
    
    _N2 = 0

    _N3 = 0

    _N4 = 0
    
    _N5 = 0
    
    _compression = 0

    _mesh = None
    
    _Nanet = 0
    _Nsnet,_Ns_n,_Ns_r = 0,0,0
    _Nsdnet = 0

    def __init__(self):
        GuidedProjectionBase.__init__(self)

        weights = {
            
        'fairness_4diff' :0,
        'fairness_diag_4diff' :0,
 
        'boundary_glide' :0, #Hui in gpbase.py doesn't work, replace here.
        'i_boundary_glide' :0,
        'boundary_z0' :0,
        'sharp_corner' : 0,
        'z0' : 0,

        'unit_edge_vec' : 0,  ## [ei, li]
        'unit_diag_edge_vec' : 0,

        'planarity' : 0,

        'orthogonal' :0,

        'Anet' : 0,  
        'Anet_diagnet' : 0,  
        
        'Snet' : 0,
        'Snet_diagnet' : 0,
        'Snet_orient' : 1,
        'Snet_constR' : 0,


        ##Note: below from geometrylab/optimization/Guidedprojection.py:
        'normal' : 0, #no use, replaced Davide's planarity way

        'edge_length' : 0,

        'area' : 0,

        'geometric' : 1, #default is True, used in constraints_equilibrium.py

        'fixed_vertices' : 1,

        'fixed_corners' : 0,

        'gliding' : 0, # Huinote: glide on boundary, used for itself boundary

        'equilibrium' : 0,
        
        'multinets_orthogonal' : 0,

        'fixed_boundary_normals': 0,        
        
        }

        self.add_weights(weights)
        
        self.switch_diagmeth = False
        
        self.is_initial = True
        
        self._glide_reference_polyline = None
        self.i_glide_bdry_crv, self.i_glide_bdry_ver = [],[]
        self.assign_coordinates = None

        self.if_uniqradius = False
        self.assigned_snet_radius = 0

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh
        self.initialization()
        
    @property
    def compression(self):
        return self._compression

    @compression.setter
    def compression(self, bool):
        if bool:
            self._compression = 2*self.get_weight('equilibrium')
            self.reinitialize = True
        elif self._compression > 0:
             self._compression = 0

    @property
    def tension(self):
        return -self._compression

    @tension.setter
    def tension(self, bool):
        if bool:
            self._compression = -2*self.get_weight('equilibrium')
            self.reinitialize = True
        elif self._compression < 0:
            self._compression = 0

    @property
    def max_weight(self):    
        return max(self.get_weight('geometric'),
                   self.get_weight('equilibrium'),
                   self.get_weight('boundary_glide'),
                   self.get_weight('planarity'),
                   
                   self.get_weight('multinets_orthogonal'),
                   
                   self.get_weight('unit_edge_vec'),
                   self.get_weight('unit_diag_edge_vec'),
                   
                   self.get_weight('orthogonal'),

                   self.get_weight('Anet'),
                   self.get_weight('Anet_diagnet'),
                   
                   self.get_weight('Snet'),
                   self.get_weight('Snet_diagnet'),

                   1)
    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self,angle):
        if angle != self._angle:
            self.mesh.angle=angle
        self._angle = angle
        
    @property
    def glide_reference_polyline(self):
        if self._glide_reference_polyline is None:
            #polylines = self.mesh.boundary_polylines() #Hui note: problem for boundary_polylines()
            polylines = self.mesh.boundary_curves(corner_split=False)[0]   
            ##print(polylines)
            
            N = 5
            try:
                for polyline in polylines:
                    polyline.refine(N)
            except:
                ###Add below for closed 1 boundary of reference mesh
                from geometrylab.geometry import Polyline
                polyline = Polyline(self.mesh.vertices[polylines],closed=True)  
                polyline.refine(N)

            self._glide_reference_polyline = polyline
        return self._glide_reference_polyline

    @glide_reference_polyline.setter
    def glide_reference_polyline(self,polyline):
        self._glide_reference_polyline = polyline        

    #--------------------------------------------------------------------------
    #                               Initialization
    #--------------------------------------------------------------------------

    def set_weights(self):
        if self.get_weight('equilibrium') != 0:
            if self.mesh.area_load != 0:
                if self.get_weight('area') == 0:
                    self.set_weight('area', 1 * self.get_weight('equilibrium'))
            else:
                self.set_weight('area', 0)
        if self.reinitialize:
            if self.get_weight('equilibrium') != 0:
                self.set_weight('edge_length', 1 * self.get_weight('equilibrium'))
                if self.mesh.area_load != 0:
                    self.set_weight('area', 1*self.get_weight('equilibrium'))
            if self.get_weight('planarity') != 0:
                self.set_weight('normal', 1*self.get_weight('planarity'))
        self.set_weight('fixed_vertices', 10 * self.max_weight)
        self.set_weight('gliding', 10 * self.max_weight)
        if self.get_weight('fixed_corners') != 0:
            self.set_weight('fixed_corners', 10 * self.max_weight)
        if self.get_weight('fixed_boundary_normals') != 0:
            self.set_weight('fixed_boundary_normals', 10 * self.max_weight)

          
    def set_dimensions(self): # Huinote: be used in guidedprojectionbase
        "X:= [Vx,Vy,Vz]"
        V = self.mesh.V
        F = self.mesh.F
        E = self.mesh.E
        N = 3*V
        N1 = N2 = N3 = N4 = N5 = N
        num_regular = self.mesh.num_regular

        Nanet = N
        Nsnet = Ns_n = Ns_r = N

        #---------------------------------------------
        if self.get_weight('planarity') != 0:
            "X += [Nx,Ny,Nz]"
            N += 3*F
            N1 = N2 = N3 = N4 = N
        if self.get_weight('equilibrium') != 0:
            "X += [edge_length, force_density, sqrt_force_density]"
            N += 3*E
            N2 = N3 = N4 = N
        if self.get_weight('area') != 0:
            "X += [Ax,Ay,Az, area]"
            N += 4*F
            N3 = N4 = N

        if self.get_weight('unit_edge_vec'): #Gnet, AGnet
            "X+=[le1,le2,le3,le4,ue1,ue2,ue3,ue4]"
            if self.get_weight('isogonal'):
                N += 16*num_regular
            else:
                "for Anet, AGnet, DGPC"
                N += 16*self.mesh.num_rrv4f4
            N5 = N
        elif self.get_weight('unit_diag_edge_vec'): #Gnet_diagnet
            "le1,le2,le3,le4,ue1,ue2,ue3,ue4 "
            N += 16*self.mesh.num_rrv4f4
            N5 = N

        if self.get_weight('Anet') or self.get_weight('Anet_diagnet'):
            N += 3*self.mesh.num_rrv4f4#3*num_regular
            Nanet = N

            
        ### Snet(_diag) project:
        if self.get_weight('Snet') or self.get_weight('Snet_diagnet'):
            num_snet = self.mesh.num_rrv4f4 
            N += 11*num_snet  
            Nsnet = N
            if self.get_weight('Snet_orient'):
                N +=4*num_snet  
                Ns_n = N
            if self.get_weight('Snet_constR'):
                N += 1
                Ns_r = N

        #---------------------------------------------
        if N1 != self._N1 or N2 != self._N2:
            self.reinitialize = True
        if N3 != self._N3 or N4 != self._N4:
            self.reinitialize = True
        if self._N2 - self._N1 == 0 and N2 - N1 > 0:
            self.mesh.reinitialize_densities()

        if N5 != self._N5:
            self.reinitialize = True
        if Nanet != self._Nanet:
            self.reinitialize = True
        if Nsnet != self._Nsnet:
            self.reinitialize = True
        if Ns_n != self._Ns_n:
            self.reinitialize = True
        if Ns_r != self._Ns_r:
            self.reinitialize = True

        #----------------------------------------------
        self._N = N
        self._N1 = N1
        self._N2 = N2
        self._N3 = N3
        self._N4 = N4
        self._N5 = N5
        self._Nanet = Nanet
        self._Nsnet,self._Ns_n,self._Ns_r = Nsnet,Ns_n,Ns_r

        self.build_added_weight() # Hui add
        
        
    def initialize_unknowns_vector(self):
        "X:= [Vx,Vy,Vz]"
        X = self.mesh.vertices.flatten('F')
        if self.get_weight('planarity') != 0:
            "X += [Nx,Ny,Nz]; len=3F"
            normals = self.mesh.face_normals()
            X = np.hstack((X, normals.flatten('F')))
        if self.get_weight('equilibrium') != 0:
            "X += [edge_length, force_density, sqrt_force_density]; len=3E"
            lengths = self.mesh.edge_lengths()
            W = self.mesh.force_densities
            X = np.hstack((X, lengths, W, np.abs(W)**0.5))
        if self.get_weight('area') != 0:
            "X += [Ax,Ay,Az, area]; len=4F"
            vector_area = self.mesh.face_vector_areas()
            face_area = np.linalg.norm(vector_area, axis=1)
            vector_area = vector_area.flatten('F')
            X = np.hstack((X, vector_area, face_area))

        if self.get_weight('unit_edge_vec'):
            _,l1,l2,l3,l4,E1,E2,E3,E4 = self.mesh.get_v4_unit_edge(rregular=True)
            X = np.r_[X,l1,l2,l3,l4]
            X = np.r_[X,E1.flatten('F'),E2.flatten('F'),E3.flatten('F'),E4.flatten('F')]

        elif self.get_weight('unit_diag_edge_vec'):
            _,l1,l2,l3,l4,E1,E2,E3,E4 = self.mesh.get_v4_diag_unit_edge()
            X = np.r_[X,l1,l2,l3,l4]
            X = np.r_[X,E1.flatten('F'),E2.flatten('F'),E3.flatten('F'),E4.flatten('F')]

        if self.get_weight('Anet') or self.get_weight('Anet_diagnet'):
            if self.get_weight('Anet'):
                if True:
                    "only r-regular vertex"
                    v = self.mesh.ver_rrv4f4
                else:
                    v = self.mesh.ver_regular
            elif self.get_weight('Anet_diagnet'):
                v = self.mesh.rr_star_corner[0]
            V4N = self.mesh.vertex_normals()[v]
            X = np.r_[X,V4N.flatten('F')]

        ### Snet-project:
        if self.get_weight('Snet') or self.get_weight('Snet_diagnet'):
            r = self.get_weight('Snet_constR')
            is_diag = False if self.get_weight('Snet') else True
            x_snet,Nv4 = self.get_snet(r,is_diag)
            X = np.r_[X,x_snet]
        #-----------------------
        
        self._X = X
        self._X0 = np.copy(X)
            
        self.build_added_weight() # Hui add

    #--------------------------------------------------------------------------
    #                       Getting (initilization + Plotting):
    #--------------------------------------------------------------------------

    def get_snet(self,is_r,is_diag=False,is_orient=True):
        """
        each vertex has one [a,b,c,d,e] for sphere equation:
            f=a(x^2+y^2+z^2)+(bx+cy+dz)+e=0
            when a=0; plane equation
            (x-x0)^2+(y-y0)^2+(z-z0)^2=R^2
            M = (x0,y0,z0)=-(b,c,d)/2a
            R^2 = (b^2+c^2+d^2-4ae)/4a^2
            unit_sphere_normal==-(2*A*Vx+B, 2*A*Vy+C, 2*A*Vz+D), 
            (note direction: from vertex to center)
            
        X += [V^2,A,B,C,D,E,a_sqrt]
        if orient:
            X += [n4,n4_sqrt], n4=-[2ax+b,2ay+c,2az+d]
        if r:
            X += [r]
        if angle
           X += [l1,l2,l3,l4,ue1,ue2,ue3,ue4]
        """
        V = self.mesh.vertices
        if is_diag:
            s0,s1,s2,s3,s4 = self.mesh.rr_star_corner
        else:
            s0,s1,s2,s3,s4 = self.mesh.rrv4f4
        S0,S1,S2,S3,S4 = V[s0],V[s1],V[s2],V[s3],V[s4]
        centers,radius,coeff,Nv4 = interpolate_sphere(S0,S1,S2,S3,S4)
        VV = np.linalg.norm(np.vstack((S0,S1,S2,S3,S4)),axis=1)**2
        A,B,C,D,E = coeff.reshape(5,-1)
        A_sqrt = np.sqrt(A)
        XA = np.r_[VV,A,B,C,D,E,A_sqrt]
        if is_orient: ##always True
            B,C,D,Nv4,n4_sqrt = self.mesh.orient(S0,A,B,C,D,Nv4)
            XA = np.r_[VV,A,B,C,D,E,A_sqrt]  
            XA = np.r_[XA, Nv4.flatten('F'),n4_sqrt]
        if is_r:
            r = np.mean(radius)
            XA = np.r_[XA,r]
        return XA, Nv4

    def get_snet_data(self,is_diag=False, ##note: combine together suitable for diagonal
                      center=False,normal=False,tangent=False,ss=False,
                      is_diag_binormal=False):
        "at star = self.rr_star"
        V = self.mesh.vertices
        if is_diag:
            s0,s1,s2,s3,s4 = self.mesh.rr_star_corner
        else:
            s0,s1,s2,s3,s4 = self.mesh.rrv4f4
            
        S0,S1,S2,S3,S4 = V[s0],V[s1],V[s2],V[s3],V[s4]
        centers,r,coeff,Nv4 = interpolate_sphere(S0,S1,S2,S3,S4)
        if self.get_weight('Snet_orient'):
            A,B,C,D,E = coeff.reshape(5,-1)
            _,_,_,Nv4,_ = self.mesh.orient(S0,A,B,C,D,Nv4)
            #Nv4 = self.mesh.get_v4_orient_unit_normal()[1][self.mesh.ind_rr_star_v4f4]
            centers = S0+r[:,None]*Nv4
        if center:
            er0 = np.abs(np.linalg.norm(S0-centers,axis=1)-r)
            er1 = np.abs(np.linalg.norm(S1-centers,axis=1)-r)
            er2 = np.abs(np.linalg.norm(S2-centers,axis=1)-r)
            er3 = np.abs(np.linalg.norm(S3-centers,axis=1)-r)
            er4 = np.abs(np.linalg.norm(S4-centers,axis=1)-r)
            #err = (er0+er1+er2+er3+er4) / 5
            err = np.sqrt(er0**2+er1**2+er2**2+er3**2+er4**2)/r
            print('radii:[min,mean,max]=','%.3f'%np.min(r),'%.3f'%np.mean(r),'%.3f'%np.max(r))
            print('Err:[min,mean,max]=','%.3g'%np.min(err),'%.3g'%np.mean(err),'%.3g'%np.max(err))
            return centers,err
        elif normal:
            # n = np.cross(C3-C1,C4-C2)
            #n = S0-centers
            n = Nv4
            un = n / np.linalg.norm(n,axis=1)[:,None]
            return S0,un
        elif is_diag_binormal:
            "only work for SSG/GSS/SSGG/GGSS-project, not for general Snet"
            n = Nv4
            un = n / np.linalg.norm(n,axis=1)[:,None]
            if is_diag:
                _,sa,sb,sc,sd = self.mesh.rrv4f4
            else:
                _,sa,sb,sc,sd = self.mesh.rr_star_corner
            t1 = (V[sa]-V[sc])/np.linalg.norm(V[sa]-V[sc], axis=1)[:,None]
            t2 = (V[sb]-V[sd])/np.linalg.norm(V[sb]-V[sd], axis=1)[:,None]
            "note works for SSG..case, since un,sa-v,sc are coplanar"
            bin1 = np.cross(un, t1)
            bin2 = np.cross(un, t2)
            bin1 = bin1 / np.linalg.norm(bin1, axis=1)[:,None]
            bin2 = bin2 / np.linalg.norm(bin2, axis=1)[:,None]
            return S0, bin1, bin2
        elif tangent:
            inn,_ = self.mesh.get_rr_vs_bounary()
            V0,V1,V2,V3,V4 = S0[inn],S1[inn],S2[inn],S3[inn],S4[inn]
            l1 = np.linalg.norm(V1-V0,axis=1)
            l2 = np.linalg.norm(V2-V0,axis=1)
            l3 = np.linalg.norm(V3-V0,axis=1)
            l4 = np.linalg.norm(V4-V0,axis=1)
            t1 = (V1-V0)*(l3**2)[:,None] - (V3-V0)*(l1**2)[:,None]
            t1 = t1 / np.linalg.norm(t1,axis=1)[:,None]
            t2 = (V2-V0)*(l4**2)[:,None] - (V4-V0)*(l2**2)[:,None]
            t2 = t2 / np.linalg.norm(t2,axis=1)[:,None]
            return V0,t1,t2
        elif ss:
            n = Nv4
            un = n / np.linalg.norm(n,axis=1)[:,None]
            return s0,un   
        return S0,S1,S2,S3,S4,centers,r,coeff,Nv4
  
    # -------------------------------------------------------------------------
    #                                 Build
    # -------------------------------------------------------------------------

    def build_iterative_constraints(self):
        self.build_added_weight() # Hui change
        
        H, r = self.mesh.iterative_constraints(**self.weights) ##NOTE: in gridshell.py
        self.add_iterative_constraint(H, r, 'mesh_iterative')
        
        H, r = self.mesh.fairness_energy(**self.weights) ##NOTE: in gridshell.py
        self.add_iterative_constraint(H, r, 'fairness')
        
        if self.get_weight('fairness_4diff'):
            pl1,pl2 = self.mesh.all_rr_continuous_polylist
            pl1.extend(pl2)
            H,r = con_fairness_4th_different_polylines(pl1,**self.weights)
            self.add_iterative_constraint(H, r, 'fairness_4diff') 
        if self.get_weight('fairness_diag_4diff'):
            pl1 = self.mesh.all_rr_diag_polylist[0][0]
            pl2 = self.mesh.all_rr_diag_polylist[1][0]
            pl1.extend(pl2)
            H,r = con_fairness_4th_different_polylines(pl1,diag=True,**self.weights)
            self.add_iterative_constraint(H, r, 'fairness_diag_4diff')     
        
        if self.get_weight('planarity'):
            #"use Hui's way not Davide's"
            H,r = con_planarity_constraints(**self.weights)  
            self.add_iterative_constraint(H, r, 'planarity')
            
        if self.get_weight('equilibrium') != 0:
            H,r = edge_length_constraints(**self.weights)
            self.add_iterative_constraint(H, r, 'edge_length')
            H,r = equilibrium_constraints(**self.weights)
            self.add_iterative_constraint(H, r,'equilibrium')
            
        if self.compression != 0 and self.get_weight('equilibrium') != 0:
            H,r = compression_constraints(self.compression,**self.weights)
            self.add_iterative_constraint(H, r, 'compression')
            
        if self.get_weight('area') != 0:
            H,r = area_constraints(**self.weights)
            self.add_iterative_constraint(H, r, 'face_vector_area')
            H,r = vector_area_constraints(**self.weights)
            self.add_iterative_constraint(H, r, 'face_area')
        
        if self.get_weight('multinets_orthogonal') !=0: 
            H,r =  con_multinets_orthogonal(**self.weights)
            self.add_iterative_constraint(H, r,'multinets_orthogonal')
            
        ###-------partially shared-used codes:---------------------------------
        if self.get_weight('unit_edge_vec'): 
            H,r = con_unit_edge(rregular=True,**self.weights)
            self.add_iterative_constraint(H, r, 'unit_edge')
        elif self.get_weight('unit_diag_edge_vec'): 
            H,r = con_unit_edge(rregular=True,**self.weights)
            self.add_iterative_constraint(H, r, 'unit_diag_edge_vec')

        if self.get_weight('boundary_z0') !=0:
            z = 0
            v = np.array([816,792,768,744,720,696,672,648,624,600,576,552,528,504,480,456,432,408,384,360,336,312,288,264,240,216,192,168,144,120,96,72,48,24,0],dtype=int)
            H,r = con_selected_vertices_glide_in_one_plane(v,2,z,**self.weights)              
            self.add_iterative_constraint(H, r, 'boundary_z0')
            
        if self.assign_coordinates is not None:
            index,Vf = self.assign_coordinates
            H,r = con_fix_vertices(index, Vf.flatten('F'),**self.weights)
            self.add_iterative_constraint(H, r, 'fix_pts')

        if self.get_weight('boundary_glide'):
            "the whole boundary"
            refPoly = self.glide_reference_polyline
            glideInd = self.mesh.boundary_curves(corner_split=False)[0] 
            w = self.get_weight('boundary_glide')
            H,r = con_alignment(w, refPoly, glideInd,**self.weights)
            self.add_iterative_constraint(H, r, 'boundary_glide')
        elif self.get_weight('i_boundary_glide'):
            "the i-th boundary"
            refPoly = self.i_glide_bdry_crv
            glideInd = self.i_glide_bdry_ver
            if len(glideInd)!=0:
                w = self.get_weight('i_boundary_glide')
                H,r = con_alignments(w, refPoly, glideInd,**self.weights)
                self.add_iterative_constraint(H, r, 'iboundary_glide')

        if self.get_weight('sharp_corner'):
            H,r = con_sharp_corner(move=0,**self.weights)
            self.add_iterative_constraint(H,r, 'sharp_corner')
            
        if self.get_weight('orthogonal'):
            H,r = con_orthogonal_midline(**self.weights)
            self.add_iterative_constraint(H, r, 'orthogonal')

        if self.get_weight('z0') !=0:
            H,r = con_glide_in_plane(2,**self.weights)
            self.add_iterative_constraint(H,r, 'z0')        
            
        if self.get_weight('Anet'):
            H,r = con_anet(rregular=True,**self.weights)
            self.add_iterative_constraint(H, r, 'Anet')
        elif self.get_weight('Anet_diagnet'):
            H,r = con_anet_diagnet(**self.weights)
            self.add_iterative_constraint(H, r, 'Anet_diagnet')

        if self.get_weight('Snet'):
            orientrn = self.mesh.new_vertex_normals()
            H,r = con_snet(orientrn,
                           is_uniqR=self.if_uniqradius,
                           assigned_r=self.assigned_snet_radius,
                           **self.weights)
            self.add_iterative_constraint(H, r, 'Snet') 
        elif self.get_weight('Snet_diagnet'):
            if 1:
                "SSG-PROJECT"
                orientrn = self.mesh.new_vertex_normals()
                H,r = con_snet(orientrn,is_diagmesh=True,
                               is_uniqR=self.if_uniqradius,
                               assigned_r=self.assigned_snet_radius,
                               **self.weights)
            else: 
                "CRPC-project"
                H,r = con_snet_diagnet(**self.weights)
            self.add_iterative_constraint(H, r, 'Snet_diag') 

        ###--------------------------------------------------------------------                
   
        self.is_initial = False   
            
        #print('-'*10)
        print(' Err_total: = ','%.3e' % np.sum(np.square(self._H*self.X-self._r)))
        #print('-'*10)
        
    def build_added_weight(self): # Hui add
        self.add_weight('mesh', self.mesh)
        self.add_weight('N', self.N)
        self.add_weight('X', self.X)
        self.add_weight('N1', self._N1)
        self.add_weight('N2', self._N2)
        self.add_weight('N3', self._N3)
        self.add_weight('N4', self._N4)
        self.add_weight('N5', self._N5)
        self.add_weight('Nanet', self._Nanet)
        self.add_weight('Nsnet', self._Nsnet)
        self.add_weight('Ns_n', self._Ns_n)
        self.add_weight('Ns_r', self._Ns_r)

    def values_from_each_iteration(self,**kwargs):
        if kwargs.get('unit_edge_vec'):
            if True:
                rr=True
            _,l1,l2,l3,l4,_,_,_,_ = self.mesh.get_v4_unit_edge(rregular=rr)
            Xi = np.r_[l1,l2,l3,l4]
            return Xi

        if kwargs.get('unit_diag_edge_vec'):
            _,l1,l2,l3,l4,_,_,_,_ = self.mesh.get_v4_diag_unit_edge()
            Xi = np.r_[l1,l2,l3,l4]
            return Xi

    def build_constant_constraints(self): #copy from guidedprojection,need to check if it works 
        self.add_weight('N', self.N)
        H, r = self.mesh.constant_constraints(**self.weights)
        self.add_constant_constraint(H, r, 'mesh_constant')
        if self.get_weight('equilibrium') > 0:
            boundary_densities_constraints(**self.weights)
        if self.get_weight('fixed_boundary_normals') > 0:
            fixed_boundary_normals_constraints(**self.weights)

    def build_constant_fairness(self): #copy from guidedprojection,need to check if it works 
        self.add_weight('N', self.N)
        K, s = self.mesh.fairness_energy(**self.weights)
        self.add_constant_fairness(K, s)
  
    def post_iteration_update(self): #copy from guidedprojection,need to check if it works 
        V = self.mesh.V
        E = self.mesh.E
        N1 = self._N1
        self.mesh.vertices[:,0] = self.X[0:V]
        self.mesh.vertices[:,1] = self.X[V:2*V]
        self.mesh.vertices[:,2] = self.X[2*V:3*V]
        if self.get_weight('equilibrium')!= 0:
            self.mesh.force_densities = self.X[N1+E:N1+2*E]
        else:
            self.mesh.force_densities = np.zeros(self.mesh.E)

    def on_reinitilize(self): #copy from guidedprojection,need to check if it works 
        self.mesh.reinitialize_force_densities()
    #--------------------------------------------------------------------------
    #                                  Results
    #--------------------------------------------------------------------------

    def vertices(self):
        V = self.mesh.V
        vertices = self.X[0:3*V]
        vertices = np.reshape(vertices, (V,3), order='F')
        return vertices

    def edge_lengths(self, initialized=False):
        if self.get_weight('equilibrium') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        E = self.mesh.E
        N1 = self._N1
        return X[N1:N1+E]

    def face_normals(self, initialized=False):
        if self.get_weight('planarity') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        V = self.mesh.V
        F = self.mesh.F
        normals = X[3*V:3*V+3*F]
        normals = np.reshape(normals, (F,3), order='F')
        return normals

    def face_vector_areas(self, initialized=False):
        if self.get_weight('area') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        F = self.mesh.F
        N2 = self._N2
        areas = X[N2:N2+3*F]
        areas = np.reshape(areas, (F,3), order='F')
        return areas

    def face_areas(self, initialized=False):
        if self.get_weight('area') == 0:
            return None
        if initialized:
            X = self._X0
        else:
            X = self.X
        F = self.mesh.F
        N2 = self._N2
        areas = X[N2+3*F:N2+4*F]
        return areas

    def force_densities(self):
        if self.get_weight('equilibrium') == 0:
            return None
        E = self.mesh.E
        N1 = self._N1
        return self.X[N1+E:N1+2*E]

    #--------------------------------------------------------------------------
    #                                Errors strings
    #--------------------------------------------------------------------------
    def make_errors(self):
        self.edge_length_error()
        self.equilibrium_error()
        self.face_areas_error()
        self.face_vector_areas_error()
        self.planarity_error() #self.face_normals_error()
        self.orthogonal_error()
        self.anet_error() 
        #self.geometric_error()

    def planarity_error(self):
        if self.get_weight('planarity') == 0:
            return None
        P = self.mesh.face_planarity()
        emean = np.mean(P)
        emax = np.max(P)
        self.add_error('planarity', emean, emax, self.get_weight('planarity'))
        print('planarity:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def orthogonal_error(self):
        if self.get_weight('orthogonal') == 0:
            return None
        _,_,t1,t2,_ = self.mesh.get_quad_midpoint_cross_vectors()
        cos = np.einsum('ij,ij->i',t1,t2)
        cos0 = np.mean(cos)
        err = np.abs(cos-cos0)
        emean = np.mean(err)
        emax = np.max(err)
        self.add_error('orthogonal', emean, emax, self.get_weight('orthogonal'))
        print('orthogonal:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def anet_error(self):
        if self.get_weight('Anet') == 0 and self.get_weight('Anet_diagnet')==0:
            return None
        if self.get_weight('Anet'):
            name = 'Anet'
            v,v1,v2,v3,v4 = self.mesh.rrv4f4
        elif self.get_weight('Anet_diagnet'):
            name = 'Anet_diagnet'
            v,v1,v2,v3,v4 = self.mesh.rr_star_corner
            
        if self.is_initial:    
            Nv = self.mesh.vertex_normals()[v]
        else:
            num = len(v)
            c_n = self._Nanet-3*num+np.arange(3*num)
            Nv = self.X[c_n].reshape(-1,3,order='F')        
        V = self.mesh.vertices
        err1 = np.abs(np.einsum('ij,ij->i',Nv,V[v1]-V[v]))
        err2 = np.abs(np.einsum('ij,ij->i',Nv,V[v2]-V[v]))
        err3 = np.abs(np.einsum('ij,ij->i',Nv,V[v3]-V[v]))
        err4 = np.abs(np.einsum('ij,ij->i',Nv,V[v4]-V[v]))
        Err = err1+err2+err3+err4
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error(name, emean, emax, self.get_weight(name))  
        print('anet:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def face_normals_error(self):
        if self.get_weight('planarity') == 0:
            return None
        N = self.face_normals()
        N0 = self.face_normals(initialized=True)
        norm = np.mean(np.linalg.norm(N, axis=1))
        Err = (np.linalg.norm(N-N0, axis=1)) / norm
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error('face_normal', emean, emax, self.get_weight('normal'))
        print('planarity:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def edge_length_error(self):
        if self.get_weight('equilibrium') == 0:
            return None
        L = self.edge_lengths()
        L0 = self.edge_lengths(initialized=True)
        norm = np.mean(L)
        Err = np.abs(L-L0) / norm
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error('edge_length', emean, emax, self.get_weight('edge_length'))
        print('edge_length:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def face_vector_areas_error(self):
        if self.get_weight('area') == 0:
            return None
        A = self.face_vector_areas()
        A0 = self.face_vector_areas(initialized=True)
        norm = np.mean(np.linalg.norm(A, axis=1))
        Err = (np.linalg.norm(A-A0, axis=1)) / norm
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error('face_vector_area', emean, emax, self.get_weight('area'))
        print('area_vector:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def face_areas_error(self):
        if self.get_weight('area') == 0:
            return None
        A = self.face_areas()
        A0 = self.face_areas(initialized=True)
        norm = np.mean(A)
        Err = (np.abs(A-A0)) / norm
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error('face_area', emean, emax, self.get_weight('area'))
        print('area:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def equilibrium_error(self):
        if self.get_weight('equilibrium') == 0:
            return None
        Err = self.mesh.equilibrium_error()
        emean = np.mean(Err)
        emax = np.max(Err)
        self.add_error('equilibrium', emean, emax, self.get_weight('equilibrium'))
        print('equilibrium:[mean,max]=','%.3g'%emean,'%.3g'%emax)

    def geometric_error(self):
        if len(self._errors) == 0:
            return None
        n = 0
        geo_mean = 0
        geo_max = 0
        if self.get_weight('planarity') != 0:
            err = self.get_error('face_normal')
            geo_mean += err[0]
            geo_max = max([geo_max, err[1]])
            n += 1
        if self.get_weight('equilibrium') != 0:
            err = self.get_error('edge_length')
            geo_mean += err[0]
            geo_max = max([geo_max, err[1]])
            n += 1
            if self.get_weight('area') != 0:
                err = self.get_error('face_vector_area')
                geo_mean += err[0]
                geo_max = max([geo_max, err[1]])
                err = self.get_error('face_area')
                geo_mean += err[0]
                geo_max = max([geo_max, err[1]])
                n += 2
        if n > 0:
            geo_mean = geo_mean / n
        self.add_error('geometric', geo_mean, geo_mean,
                       self.get_weight('geometric'))

    def geometric_error_string(self):
        return self.error_string('geometric')

    def equilibrium_error_string(self):
        return self.error_string('equilibrium')

    def planarity_error_string(self):
        return self.error_string('planarity')

    def orthogonal_error_string(self):
        return self.error_string('orthogonal')
    
    def anet_error_string(self):
        return self.error_string('Anet')

    #--------------------------------------------------------------------------
    #                                   Utilities
    #--------------------------------------------------------------------------

    def axial_forces(self):
        if self.equilibrium != 0:
            return self.mesh.axial_forces()
        else:
            return np.zeros(self.mesh.E)

    def force_resultants(self):
        return self.mesh.force_resultants()

    def applied_loads(self):
        return self.mesh.applied_loads()
