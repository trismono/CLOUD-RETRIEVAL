#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np
from numpy.linalg import inv

def l2m_retrieval(rad_meas,rad_fwd,sy,sa,kmat,xii,xa,gamma):
    # Identity matrix (n state) for scaling
    D = np.identity(len(sa))

    # Levenberg-Marquardt (Rodgers 2000)
    # xi+1 = xi + [(Sa-1 + Kit Sy-1 Ki + gamma D-1]-1 {Kit Sy-1 [y - F(xi)] - Sa-1 [xi - xa]}
    x = xii + np.matmul(inv(inv(sa) + np.matmul(kmat.T,np.matmul(inv(sy),kmat)) + (gamma*inv(D))),(np.matmul(kmat.T,np.matmul(inv(sy),rad_meas-rad_fwd)) - np.matmul(inv(sa),xii-xa)))

    # Covariance matrix of the solution
    # Sx = (Kt Sy-1 K + Sa-1)-1
    sx = inv(np.matmul(kmat.T,np.matmul(inv(sy),kmat)) + inv(sa))
    
    # Gain-matrix
    # G = (Kt Sy-1 K + Sa-1)-1 Kt Sy-1 = Sx Kt Sy-1
    G = np.matmul(sx,np.matmul(kmat.T,inv(sy)))

    # Averaging kernel
    # A = G K
    A = np.matmul(G,kmat)   
    
    # Cost function
    # chi2 = [y − F(x)]t Sy−1 [y − F(x)] + [xi − xa]T Sa−1 [xi − xa]
    chi2 = np.matmul((rad_meas-rad_fwd).T,np.matmul(inv(sy),rad_meas-rad_fwd)) + np.matmul((x-xa).T,np.matmul(inv(sa),x-xa)) 
    
    # return output
    return(x,sx,G,A,chi2)