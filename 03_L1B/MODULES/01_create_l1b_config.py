#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

# import standard modules
import os

# limit openblas thread nuumber
os.environ["OPENBLAS_NUM_THREADS"] = "3"

import numpy as np
from configparser import ConfigParser

# get and define working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define path
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"

# +++++++++++++++++++
# define SNR of SMART
# +++++++++++++++++++
# here SMART SNR is derived from modis
# it's assumed that SMART noises are two times higher than MODIS
# in the future, real SNR should be used

# define central modis wavelength
lambda_modis = np.array([620 + 670	,
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
snr_modis = np.array([128,
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

# define SMART SNR
snr_smart = snr_modis/2

# convert to list and to string
lambda_modis = ', '.join((lambda_modis.astype(str)).tolist())
snr_smart = ', '.join((snr_smart.astype(str)).tolist())

# ++++++++++++++++++++++
# define isrf parameters
# ++++++++++++++++++++++
# SMART ISRF (slit function) parameters
fwhm = 3                # full with at half maximum [nm]
nisrf_grid = 201        # number of ISRF grids
lambda_offset = 5       # offset from the center [nm]

# define filename
filename = input_dir + "l1b_config.dat"

# create input directory
try:
    os.makedirs(input_dir)
except OSError:
    pass

# create config file
config = ConfigParser()

# libradtran
config["SNR"] = {
    "wavelength": "%s" %lambda_modis,
    "snr_smart": "%s" %snr_smart
    }

# number of core
config["ISRF"] = {
    "fwhm": "%.2f" %fwhm,
    "nisrf_grid": "%i" %nisrf_grid,
    "lambda_offset": "%.2f" %lambda_offset,
    }

# +++++++++++++++++++
# writing config file
# +++++++++++++++++++
with open(filename, "w") as f:
    config.write(f)