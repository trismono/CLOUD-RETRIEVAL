#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np

def diagnostic(stat_file,pix_id,message,niter,chi2,x,gain_file,lambda_meas,G):    
    # open statistic file
    with open(stat_file, 'w+') as fid:
        # write a header :: variable names and units
        fid.write("  pix_id    conv_id     n_iter         chi2        tau   reff (\mu) \n")
        
        # write retrieval outputs
        fid.write("%8i %10i %10i %12.3f %10.3f %12.3f" %(pix_id,       # pixel idntifier
                                                         message,      # convergence diagnostic
                                                         niter,        # number of iteration
                                                         chi2,         # cost function
                                                         x[0],         # tau
                                                         x[1]))        # reff
    
    # create array for gain :: lambda, gain_tau, gain_reff
    data = np.zeros(shape=(len(lambda_meas),3), dtype = float)
    
    # assign data to array
    data[:,0] = lambda_meas     # lambda
    data[:,1:] = G.T            # gain :: original dimension = nstate * n_meas
    
    # write gain matrix into file :: lambda, gain_tau, gain_reff
    np.savetxt(gain_file, data, fmt=["%10.3f", "%20e", "%20e"])    
    
    
    
    
    