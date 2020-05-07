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

import numpy as np
from config_module import config_module
import itertools

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define dirs
data_dir = os.path.realpath(base_dir + "../../DATA") + "/" 
out_dir = os.path.realpath(base_dir + "../OUTPUT") + "/" 

# create dir
try:
    os.makedirs(out_dir)
except OSError:
    pass

# +++++++++++++++++
# measurement setup
# +++++++++++++++++
sza = 30.0
phi0 = 100.0
phi = 0
zout = 3
doy = 85    # day of year
albedo_file = data_dir + "albedo_sea_water.dat"
ztop = 2.0
zbase = 1.5
tau550 = np.array([4,6,8])
reff = np.array([8,10,12]) # micron
aerosol_haze = "4"
aerosol_season = "1"
lambda0 = 400
lambda1 = 2000

# ++++++++++++++++++
# create combination
# ++++++++++++++++++
comb = itertools.product(tau550,reff)
comb = np.array(list(comb))

# ++++++++++++++++++
# create cases
# ++++++++++++++++++
# define dimension
nobs = len(comb)

# loop overr pixels
for i in range(nobs):
    # define filename
    outfile = out_dir + "meas_setup_%i.dat" %(i+1)
    
    # call config module
    config_module(filename=outfile,sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,\
                  albedo_file=albedo_file,ztop=ztop,zbase=zbase,tau550=comb[i,0],\
                  reff=comb[i,1],aerosol_season=aerosol_season,aerosol_haze=aerosol_haze,\
                  lambda0=lambda0,lambda1=lambda1)
