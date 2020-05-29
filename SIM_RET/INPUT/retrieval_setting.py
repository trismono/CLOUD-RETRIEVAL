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

# limit the number of cores
os.environ["OPENBLAS_NUM_THREADS"] = "3"      # openblas thread number (problem for Python 3)

from configparser import ConfigParser

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define data dir
data_dir = os.path.realpath(base_dir + "../../DATA") + "/"

# define outfile
filename = base_dir + "retrieval_setting.ini"

# define retrieval parameters
tau0 = 8
reff0 = 10
weighting_tau0 = tau0*2
weighting_reff = reff0*2
ztop = 2.5
zbase = 2.0
aerosol_haze = "4"
aerosol_season = "1"
lambda0_fwd = 550
lambda1_fwd = 1700
lambda0_ret_vnir = 630
lambda1_ret_vnir = 670
lambda0_ret_swir = 1500
lambda1_ret_swir = 1660
albedo_file = data_dir + "albedo_sea_water.dat"
iter_max = 15
noise = 1           # 0 = without noise, 1 = with noise
rad_noise = 0.03

# create config file
config = ConfigParser()

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
