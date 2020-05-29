#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

from configparser import ConfigParser

def config_module(filename,sza,phi0,phi,zout,doy,albedo_file,ztop,zbase,\
                  reff,tau550,aerosol_season,aerosol_haze,lambda0,lambda1):
    
    # create config file
    config = ConfigParser()
    
    # sun and measurement geometry
    config["GEOMETRY"] = {
            "sza": "%.2f" %sza,
            "phi0": "%.2f" %phi0,
            "phi": "%.2f" %phi,
            "zout": "%.2f" %zout,
            "doy": "%i" %doy
            }
    
    # albedo file
    config["ALBEDO_FILE"] = {
            "albedo_file": albedo_file
            }
    
    # cloud profile (simplified cloud model)
    config["CLOUD_PROFILE"] = {
            "ztop": "%.2f" %ztop,
            "zbase": "%.2f" %zbase,
            "reff": "%.2f" %reff,
            "tau550": "%.2f" %tau550
            }
    
    # aerosol
    config["AEROSOL"] = {
            "aerosol_season": aerosol_season,
            "aerosol_haze": aerosol_haze
            }

    # wavelengths
    config["WAVELENGTH"] = {
            "lambda0": "%.2f" %lambda0,
            "lambda1": "%.2f" %lambda1
            }
    
    # writing file
    with open(filename, "w") as f:
        config.write(f)
