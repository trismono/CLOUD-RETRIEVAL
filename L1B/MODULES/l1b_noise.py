#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np

def l1b_noise(wavelength, spectra):
    # define spectra dimension
    nspectra = len(spectra)
    
    # +++++++++++++++++++
    # take snr from modis
    # +++++++++++++++++++
    # define modis wavelength
    lambda_mod = np.array([620 + 670	,
                           841 + 876	,
                           459 + 479	,
                           545 + 565	,
                           1230 + 1250,
                           1628 + 1652,
                           2105 + 2155,
                           405 + 420	,
                           438 + 448	,
                           438 + 493	,
                           526 + 536	,
                           546 + 556	,
                           662 + 672	,
                           673 + 683	,
                           743 + 753	,
                           862 + 877	,
                           890 + 920	,
                           931 + 941	,
                           915 + 965	])/2
    
    # define MODIS SNR
    snr_mod = np.array([128,
                        201,
                        243,
                        228,
                        74,
                        274,
                        110,
                        880,
                        838,
                        802,
                        754,
                        750,
                        910,
                        1087,
                        586,
                        516,
                        167,
                        57,
                        250])
    
    # define SMART SNR : half of MODIS SNR
    snr_smart = snr_mod/2
    
    # +++++++++++++
    # fix the order
    # +++++++++++++
    # find index
    idx = np.argsort(lambda_mod)
    
    # reorder lambda and snr
    lambda_mod = lambda_mod[idx]
    snr_mod = snr_mod[idx]          # it's not used thereafter
    snr_smart = snr_smart[idx]

    # +++++++++++++    
    # interpolation
    # +++++++++++++
    # define snr at SMART wavelengths
    snr_smart_full = np.interp(wavelength,lambda_mod,snr_smart)
    
    # define SMART noise : noise = signal/snr
    noise_smart = spectra/snr_smart_full

    # +++++++++++++++++++++++
    # random number generator 
    # +++++++++++++++++++++++
    # define and and std
    mean = 0
    std = 1
    
    # generate random number    
    seed = np.random.normal(mean,std,nspectra)
    
    # ++++++++++++++++++++
    # generate SMART noise
    # ++++++++++++++++++++
    # define error : error = noise*seed
    error = noise_smart*seed
    
    # define noisy spectra
    noisy_spectra = spectra+error
    
    # return outputs
    return(noisy_spectra, error)
    
    



