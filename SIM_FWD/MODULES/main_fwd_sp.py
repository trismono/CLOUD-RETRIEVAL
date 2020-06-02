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
os.environ["OPENBLAS_NUM_THREADS"] = "3"      # limit number of thread (problem for numpy - Python 3)
import numpy as np
import glob
from fwd_wrapper import fwd_wrapper

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define paths
case_dir = os.path.realpath(base_dir + "../../CREATE_CASE/OUTPUT") + "/"
inout_dir = os.path.realpath(base_dir + "../INOUT") + "/"
cloud_dir = os.path.realpath(base_dir + "../CLOUD") + "/"

# create input-output directory
try:
    os.makedirs(inout_dir)
    os.makedirs(cloud_dir)
except OSError:
    pass

# +++++++++++++++++++++
# running forward model
# +++++++++++++++++++++
# define number of pixels
config_list = glob.glob(case_dir + "meas_setup_*.dat")
config_list = np.sort(config_list)
nfiles = len(config_list) 

# loop over pixels
for i in range(nfiles):
    # define pixel id
    pix_id = i+1    
    
    # run forward model
    fwd_wrapper(base_dir=base_dir, pix_id=pix_id)
