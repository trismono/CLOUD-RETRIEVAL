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
os.environ["OPENBLAS_NUM_THREADS"] = "3"      # openblas thread number (problem for Python 3)
import numpy as np
import glob
import concurrent.futures

# import local routines
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

# ++++++++++++++++++++++
# multiprocessing module
# ++++++++++++++++++++++
# define number of process ::  change this number accordingly
# when running the code in the local machine, optimum performance
# is normally n_proc = number of CPU (core) - 1
# when running the code on the server, please beware of 
# processes by others and memory consumption
n_proc = 3

# define number of pixels
config_list = glob.glob(case_dir + "meas_setup_*.dat")
nfiles = len(config_list) 

# define iterable arguments for multiprocessing module
base_dir_iter = np.repeat(base_dir,nfiles)
pix_id_iter = np.arange(1,nfiles+1,1)

# initialize multiprocessing module
with concurrent.futures.ProcessPoolExecutor(max_workers=n_proc) as executor:
    
    # running retrieval l2m in parallel
    executor.map(fwd_wrapper, base_dir_iter, pix_id_iter) 
