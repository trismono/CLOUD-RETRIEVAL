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
from configparser import ConfigParser

def l1b_convolution(config_file, input_file, spectral_grid, output_file):
    # ++++++++++++++++++
    # read input spectra
    # ++++++++++++++++++
    # fwd output
    data_fwd = np.loadtxt(input_file, skiprows=1)
    
    # define dimension
    ntype = data_fwd.shape[1]
    
    # define lambda min and max
    lambda_min, lambda_max = np.min(data_fwd[:,0]), np.max(data_fwd[:,0])

    # +++++++++++    
    # create ISRF
    # +++++++++++
    # open config file
    config = ConfigParser()
    config.read(config_file)
    
    # parse config
    fwhm = float(config['ISRF']['fwhm'])
    nisrf_grid = int(config['ISRF']['nisrf_grid']) 
    lambda_offset = float(config['ISRF']['lambda_offset']) 
    
    # define static parameters
    std = fwhm/2.355    # based on the relatinship between FWHM and stdev
    dlambda =  np.linspace(-lambda_offset,lambda_offset,nisrf_grid)     # ISRF wavelength grid

    # create isrf and normalization   
    isrf = norm.pdf(dlambda,0,std)
    isrf = isrf/np.trapz(isrf,dlambda)
    
    # ++++++++++++++++++++++++++++++++++
    # trim spectral grid for convolution
    # ++++++++++++++++++++++++++++++++++
    # find index from lower boundary
    idx = np.where(spectral_grid <= lambda_min+lambda_offset)
    
    # trim lower boundary
    spectral_grid = np.delete(spectral_grid, idx)
    
    # find index from upper boundary
    idx = np.where(spectral_grid >= lambda_max-lambda_offset)
    
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
    
    # ++++++++++++
    # write output
    # ++++++++++++
    # open output file
    file = open(output_file, "w")
    file.write(" lambda   Irr_dw       Irr_up       Rad_dw       Rad_up       Act_dw       Act_up       Irr_dw_dir   Irr_up_diff  Act_dw_dir   Act_up_diff\n")

    # define the format
    output_format = "%8.3f  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n"
    
    # write output
    for i in range(nspectra):
        file.write(output_format %(spectra_data[i,0], spectra_data[i,1], spectra_data[i,2], \
                                   spectra_data[i,3], spectra_data[i,4], spectra_data[i,5], \
                                   spectra_data[i,6], spectra_data[i,7], spectra_data[i,8], \
                                   spectra_data[i,9], spectra_data[i,10]))
    
    # close file
    file.close()   
