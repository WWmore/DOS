# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 19:01:26 2019

@author: hwang
"""
__author__ = 'Hui'
#------------------------------------------------------------------------------
import numpy as np
#------------------------------------------------------------------------------


class Basis(object):
    
    def __init__(self, vectors, **kwargs):
        
        self.N = vectors

        self.V0 = kwargs.get('anchors',  None) # circles center

        self.ri = kwargs.get('r', 1)
        
        
    def orthogonalPlane3Vectors(self):
        "get three unit vectors E1,E2,E3 on planes that are orthogonal to vectors N"
        if self.V0 is not None: 
            P1, P2, P3 = [], [], []
            V0 = self.V0
        else:
            E1, E2 ,E3 = [], [], []
        N = self.N
        num = len(N)
        
        arrR = self.ri * np.ones(num)
        
        for i in range(num):
            #x,y,z = round(N[i,0],6),round(N[i,1],6),round(N[i,2],6)
            x,y,z = N[i,0],N[i,1],N[i,2]            
            if x==0 or y==0 or z==0:
                e1, e2, e3 = [0,0,0],[0,0,0],[0,0,0]
                j = 0 if x==0 else 1 if y==0 else 2
                e1[j],e2[j],e3[j]=1,2,3
                e1[j-1],e1[j-2]=-N[i,j-2],N[i,j-1]
                e2[j-1],e2[j-2]=-N[i,j-2],N[i,j-1]
                e3[j-1],e3[j-2]=-N[i,j-2],N[i,j-1]
                if (y==0 and z==0) or (x==0 and z==0) or (x==0 and y==0):  
                    e1, e2, e3 = [0,0,0],[0,0,0],[0,0,0]    
                    j = 0 if (y==0 and z==0) else 1 if  (x==0 and z==0) else 2
                    e1[j-1]=1
                    e2[j-2]=1
                    e3[j-1],e3[j-2]=-1,-1
            else:       
                e1 = [0,-z,y]
                e2 = [-z,0,x]
                e3 = [-y,x,0]
            e1,e2,e3 = np.array(e1),np.array(e2),np.array(e3)
            eps = np.finfo(float).eps
            e1=e1/(np.linalg.norm(e1)+eps)
            e2=e2/(np.linalg.norm(e2)+eps)
            e3=e3/(np.linalg.norm(e3)+eps)
            if self.V0 is not None: 
                P1.append(V0[i]+e1*arrR[i])
                P2.append(V0[i]+e2*arrR[i])
                P3.append(V0[i]+e3*arrR[i])
            else:
                E1.append(e1)
                E2.append(e2)
                E3.append(e3)
                
        if self.V0 is not None: 
            return P1,P2,P3
        else:        
            return np.array(E1),np.array(E2),np.array(E3) 
    
    def orthogonalBasis(self):
        "get another two unit vectors together with it to form a Frame"
        E1,_,_ = self.orthogonalPlane3Vectors()
        eps = np.finfo(float).eps
        E0 = self.N / (np.linalg.norm(self.N, axis=1, keepdims=True)+eps)
        E2 = np.cross(E0, E1)
        E2 = E2 / (np.linalg.norm(E2, axis=1, keepdims=True)+eps)
        return E0, E1, E2

            
