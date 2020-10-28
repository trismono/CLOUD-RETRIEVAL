#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

# import modules
import os

# limit openblas thread number
os.environ["OPENBLAS_NUM_THREADS"] = "3"

from configparser import ConfigParser

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define path
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"
data_dir = os.path.realpath(base_dir + "../../DATA") + "/"

# define outfile
filename = input_dir + "error_config.dat"

# define retrieval parameters
tau0 = 8                                                # initial guess - optical thickness
reff0 = 10                                              # initial guess - effective radius
weighting_tau0 = tau0*2                                 # covariance matrix - optical thickness
weighting_reff = reff0*2                                # covariance matrix - effective radius
ztop = 2.5                                              # cloud top altitude
zbase = 2.0                                             # cloud base altitude
aerosol_haze = "4"                                      # aerosol haze - see libradtran user guideline
aerosol_season = "1"                                    # aerosol season - see libradtran user guideline
lambda0_fwd = 550                                       # start of simulated wavelengths
lambda1_fwd = 1700                                      # end of simulated wavelengths
lambda0_ret_vnir = 630                                  # start of vnir fitting window
lambda1_ret_vnir = 670                                  # end of vnir fitting window
lambda0_ret_swir = 1500                                 # start of swir fitting window
lambda1_ret_swir = 1660                                 # end of swir fitting window
albedo_file = data_dir + "albedo_sea_water.dat"         # albedo file
iter_max = 15                                           # number of maximum iteration
noise = 1                                               # 0 = without noise, 1 = with noise
rad_noise = 0.03                                        # 0.03 = 3%
libradtran_dir = "/deos/trismonock/libRadtran-2.0.2/"   # this must be adjusted locally
n_core = 3                                              # number of core used 

# create config file
config = ConfigParser()

# libradtran directory
config["LIBRADTRAN"] = {
        "LIBRADTRAN_DIR": "%s" %libradtran_dir,
        }

# number of core
config["CORE"] = {
        "n_core": "%i" %n_core,
        }

# initial guess
config["FLAG"] = {
        "noise": "%i" %noise,
        }

# initial guess
config["PRIOR_VALUE"] = {
        "tau0": "%10e" %tau0,
        "reff0": "%10e" %reff0
        }

# wighting
config["WEIGHTING"] = {
        "tau": "%10e" %weighting_tau0,
        "reff": "%10e" %weighting_reff
        }

# cloud
config["CLOUD_PROFILE"] = {
        "ztop": "%10e" %ztop,
        "zbase": "%10e" %zbase
        }

# aerosol
config["AEROSOL"] = {
        "aerosol_season": aerosol_season,
        "aerosol_haze": aerosol_haze
        }

# wavelength fwd sim
config["WAVELENGTH_FWD"] = {
        "lambda0": "%.2f" %lambda0_fwd,
        "lambda1": "%.2f" %lambda1_fwd
        }

# wavelength ret
config["WAVELENGTH_RET"] = {
        "lambda0_ret_vnir": "%10e" %lambda0_ret_vnir,
        "lambda1_ret_vnir": "%10e" %lambda1_ret_vnir,
        "lambda0_ret_swir": "%10e" %lambda0_ret_swir,
        "lambda1_ret_swir": "%10e" %lambda1_ret_swir
        }

# albedo
config["ALBEDO"] = {
        "albedo_file": albedo_file
        }

# number of iteration
config["ITERATION"] = {
        "iter_max": "%i" %iter_max
        }

# radiance noise
config["NOISE"] = {
        "rad_noise": "%10e" %rad_noise
        }

# writing file
with open(filename, "w") as f:
    config.write(f)
