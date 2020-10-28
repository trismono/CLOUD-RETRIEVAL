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

def l2m_retrieval(rad_meas, rad_fwd, sy, sa, kmat, xi, xa, gamma):
    
    """
    retrievals are carried out using bayesian optimal estimation
    here we use a modified Gauss-Newton Levenberg-Marquardt method
    for more detail please refer to Rodgers (2000)
    """
    
    # Gauss-Newton iterative method
    # xi+1 = xa + [gamma Sa-1 + KiT Sy-1 Ki]-1 KiT Sy-1 [y - F(x)i + Ki (xi - xa)]
    x = xa + np.matmul(inv(gamma*inv(sa) + np.matmul(kmat.T,np.matmul(inv(sy),kmat))),np.matmul(kmat.T,np.matmul(inv(sy),(rad_meas-rad_fwd + np.matmul(kmat,xi-xa)))))

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
    
    # return output:
    # x     = state vector
    # sx    = posterior covariance matrix
    # G     = gain matrix
    # A     = averaging kernel
    # chi2  = cost function
    return(x, sx, G, A, chi2)