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
os.environ["OPENBLAS_NUM_THREADS"] = "3"      # limit openblas thread number (problem for numpy - Python 3)
import numpy as np
import glob
from l2m_wrapper import l2m_wrapper

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# +++++++++++++++++++++   
# define and create dir
# +++++++++++++++++++++
# define dirs
case_dir = os.path.realpath(base_dir + "../../CREATE_CASE/OUTPUT") + "/"    # to get measurement geometry
dummy_dir = os.path.realpath(base_dir + "../DUMMY") + "/"       # to store inp-out fwd sim and spectra
out_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"        # to store output

# create dirs
try:
    os.makedirs(dummy_dir)
    os.makedirs(out_dir)
except OSError:
    pass

# ++++++++++++++++++++++++
# serial processing module
# ++++++++++++++++++++++++
# define number of pixels
config_list = glob.glob(case_dir + "meas_setup_*.dat")
config_list = np.sort(config_list)
nfiles = len(config_list)

# loop over pixels
#for i in range(nfiles):
for i in [1]:
    # define pix_id
    pix_id = i+1
    
    # call l2m wrapper
    l2m_wrapper(base_dir=base_dir,pix_id=pix_id)
