#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np

def kmat_module(data_ref,data_ptau,data_preff,tau1,reff1,tau2,reff2,lambda_meas,rad_meas,sy,outfile):
    
    # ++++++++++
    # parse data
    # ++++++++++
    # 0 = lambda, 1 = irradiance down, 2 = radiance up
    # reference spectra
    rad_ref = data_ref[:,2]
    
    # perturbed tau
    rad_ptau = data_ptau[:,2]
    
    # perturbed reff
    rad_preff = data_preff[:,2]
    
    # define dimension
    nspectra = len(rad_ref)
     
    # +++++++++++++++++++++++++++++
    # calculate jacobian derivative
    # +++++++++++++++++++++++++++++
    # tau
    kmat_tau = (rad_ref-rad_ptau)/(tau1-tau2)
    
    # reff
    kmat_reff = (rad_ref-rad_preff)/(reff1-reff2)

    # ++++++++++++
    # combine data
    # ++++++++++++
    # define dimension column for wavelength, spectra_meas, sy, spectra_mod, kmat_tau, kmat_reff
    ncol = 6    
    
    # create array
    kmat_array =  np.zeros(shape=(nspectra,ncol)) 
    
    # assign data to array
    kmat_array[:,0] = lambda_meas
    kmat_array[:,1] = rad_meas
    kmat_array[:,2] = sy
    kmat_array[:,3] = rad_ref
    kmat_array[:,4] = kmat_tau
    kmat_array[:,5] = kmat_reff

    # ++++++++++++
    # write output
    # ++++++++++++    
    # open file
    file = open(outfile, "w")
    
    # write header
    file.write("%10s %20s %20s %20s %20s %20s\n" %("lambda", "Rad_up_meas", "sy", "Rad_up_mod", "Kmat_tau", "Kmat_reff"))

    # define format
    format_str = "%10.3f %20e %20e %20e %20e %20e\n"
    
    # loop over spectra
    for i in range(nspectra):
        # write output
        file.write(format_str %(kmat_array[i,0], kmat_array[i,1], kmat_array[i,2], kmat_array[i,3], kmat_array[i,4], kmat_array[i,5]))
    
    # close file
    file.close()     
    
    # return output
    return(kmat_array)