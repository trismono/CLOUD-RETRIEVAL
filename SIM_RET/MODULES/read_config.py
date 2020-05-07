#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

from configparser import ConfigParser

def read_ret_config(config_file):
    # open and parse config file
    config = ConfigParser()
    config.read(config_file)
    
    # albedo    
    albedo_file = config['ALBEDO']['albedo_file']

    # cloud profile    
    ztop = config['CLOUD_PROFILE']['ztop']
    zbase = config['CLOUD_PROFILE']['zbase']
    
    # aerosol
    aerosol_season = config['AEROSOL']['aerosol_season']
    aerosol_haze = config['AEROSOL']['aerosol_haze']
    
    # wavelength fwd
    lambda0_fwd = config['WAVELENGTH_FWD']['lambda0']
    lambda1_fwd = config['WAVELENGTH_FWD']['lambda1']

    # wavelength ret
    lambda0_ret_vnir = config['WAVELENGTH_RET']['lambda0_ret_vnir']
    lambda1_ret_vnir = config['WAVELENGTH_RET']['lambda1_ret_vnir']
    lambda0_ret_swir = config['WAVELENGTH_RET']['lambda0_ret_swir']
    lambda1_ret_swir = config['WAVELENGTH_RET']['lambda1_ret_swir']

    # prior value
    tau0 = config['PRIOR_VALUE']['tau0']
    reff0 = config['PRIOR_VALUE']['reff0']
    
    # weighting
    weighting_tau0 = config['WEIGHTING']['tau']
    weighting_reff0 = config['WEIGHTING']['reff']

    # iteration
    iter_max = config['ITERATION']['iter_max']
    
    # flag
    flag_noise = config['FLAG']['noise']

    # noise
    rad_noise = config['NOISE']['rad_noise']
    
    # return outputs
    return(albedo_file,ztop,zbase,tau0,reff0,weighting_tau0,weighting_reff0,\
           aerosol_season,aerosol_haze,lambda0_fwd,lambda1_fwd,lambda0_ret_vnir,\
           lambda1_ret_vnir,lambda0_ret_swir,lambda1_ret_swir,flag_noise,rad_noise,iter_max)

def read_geometry_config(config_file):
    # open and parse config file
    config = ConfigParser()
    config.read(config_file)
    
    # geometry
    sza = config['GEOMETRY']['sza']
    phi0 = config['GEOMETRY']['phi0']
    phi = config['GEOMETRY']['phi']
    zout = config['GEOMETRY']['zout']
    doy = config['GEOMETRY']['doy']
    
    # return outputs
    return(sza,phi0,phi,zout,doy)
    