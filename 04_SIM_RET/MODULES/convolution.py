#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np
from scipy.stats import norm 

def convolution(data_fwd, spectral_grid):
    # ++++++++++++++++++
    # read input spectra
    # ++++++++++++++++++    
    # define dimension
    ntype = data_fwd.shape[1]
    
    # define lambda min and max
    lambda_min, lambda_max = np.min(data_fwd[:,0]), np.max(data_fwd[:,0])

    # +++++++++++    
    # create ISRF
    # +++++++++++
    # define static parameters
    fwhm = 3            # as written in my thesis :D
    std = fwhm/2.355    # based on the relatinship between FWHM and stdev (source: Wikipedia)
    nlambda = 201       # fine enough
    dlambda0 = -5       # offset lower boundary
    dlambda1 = 5        # offset upper boundary
    dlambda =  np.linspace(dlambda0,dlambda1,nlambda)     # ISRF wavelength grid

    # create isrf (+renormalization)    
    isrf = norm.pdf(dlambda,0,std)
    isrf = isrf/np.trapz(isrf,dlambda)
    
    # ++++++++++++++++++++++++++++++++++
    # trim spectral grid for convolution
    # ++++++++++++++++++++++++++++++++++
    # find index from lower boundary
    idx = np.where(spectral_grid <= lambda_min)
    
    # trim lower boundary
    spectral_grid = np.delete(spectral_grid, idx)
    
    # find index from upper boundary
    idx = np.where(spectral_grid >= lambda_max)
    
    # trim upper boundary
    spectral_grid = np.delete(spectral_grid, idx)    

    # ++++++++++++++++++++
    # spectral convolution
    # ++++++++++++++++++++        
    # define dimension
    nspectra = len(spectral_grid)

    # create empty array
    spectra_data = np.zeros(shape=(nspectra,ntype))
    
    # lopp over spectra
    for i in range(nspectra):
        # define grid to be interpolated
        xdum = spectral_grid[i]+dlambda
        
        # assign wavelength
        spectra_data[i,0] = spectral_grid[i]

        # loop over quantities        
        for j in range(1,ntype):
            # interpolation to isrf wavelength grid
            dummy = np.interp(xdum,data_fwd[:,0],data_fwd[:,j])
    
            # convolution  
            spectra_data[i,j] = np.trapz(dummy*isrf,xdum)/np.trapz(isrf,xdum)
    
    # return output
    return(spectra_data)
