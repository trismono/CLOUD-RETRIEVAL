#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np
from configparser import ConfigParser

def l1b_noise(config_file, wavelength, spectra):
    # define spectra dimension
    nspectra = len(spectra)
    
    # ++++++++++++++++++++
    # read snr information
    # ++++++++++++++++++++    
    # open config file
    config = ConfigParser()
    config.read(config_file)
    
    # parse config
    snr_wavelength = np.array((config['SNR']['wavelength']).split(","), dtype=float)
    snr_smart = np.array((config['SNR']['snr_smart']).split(","), dtype=float)
    
    # +++++++++++++
    # fix the order
    # +++++++++++++
    # find index
    idx = np.argsort(snr_wavelength)
    
    # reorder lambda and snr
    snr_wavelength = snr_wavelength[idx]
    snr_smart = snr_smart[idx]

    # +++++++++++++    
    # interpolation
    # +++++++++++++
    # define snr at SMART wavelengths
    snr_smart_full = np.interp(wavelength,snr_wavelength,snr_smart)
    
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
    
    



