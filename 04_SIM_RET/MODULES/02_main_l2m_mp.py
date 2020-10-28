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

# limit openblas thread number
os.environ["OPENBLAS_NUM_THREADS"] = "3"

import numpy as np
import glob
import concurrent.futures
from configparser import ConfigParser

from l2m_wrapper import l2m_wrapper

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define path
case_dir = os.path.realpath(base_dir + "../../01_CREATE_CASE/OUTPUT") + "/" # measurement geometry
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"                   # to store inp-out uvspec and cloud file
out_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"                    # to store output

# create dirs
try:
    os.makedirs(input_dir)
    os.makedirs(out_dir)
except OSError:
    pass

# print statement
dashed = 60 * "+"
print("Info         |", dashed)
print("Info         | End-to-End Simulator Level-2 Retrieval Module")
print("Info         | Author       : Trismono C. Krisna")
print("Info         | Institution  : SRON Netherlands Institute for Space Research")
print("Info         | Version      : v2.2")
print("Info         | Date         : 27 October 2020")
print("Info         |", dashed)

# ++++++++++++
# read setting
# ++++++++++++
# define config file
config_file = input_dir + "l2m_config.dat"

# open config file
config = ConfigParser()
config.read(config_file)

# parse config
libradtran_dir = config['LIBRADTRAN']['LIBRADTRAN_DIR']
noise = config['FLAG']['noise']
n_core = int(config['CORE']['n_core'])
tau0 = float(config['PRIOR_VALUE']['tau0'])
reff0 = float(config['PRIOR_VALUE']['reff0'])
weighting_tau = float(config['WEIGHTING']['tau'])
weighting_reff = float(config['WEIGHTING']['reff'])
iter_max = int(config['ITERATION']['iter_max'])
rad_noise = float(config['NOISE']['rad_noise'])*100

# convert noise to string
if noise == "0":
    noise_str = "without noise"
elif noise == "1":
    noise_str = "with noise"

# print statement
print("Info         | ")
print("Info         | Noise realisation :", noise_str)
print("Info         | Initial guess of tau :", tau0)
print("Info         | Initial guess of reff :", reff0, "micron")
print("Info         | Variance of tau :", weighting_tau)
print("Info         | Variance of reff :", weighting_reff, "micron")
print("Info         | Number of maximum iteration :", iter_max)
print("Info         | Radiance noise :", rad_noise, "%")
print("Info         | Number of core :", n_core)
print("Info         | LibRadtran dir :", libradtran_dir)
print("Info         | ") 

# ++++++++++++++++++++++
# multiprocessing module
# ++++++++++++++++++++++
# define number of process ::  change this number accordingly
# when running the code in the local machine, optimum performance
# is normally n_proc = number of CPU (core) - 1
# when running the code on the server, please beware of 
# processes by others and memory consumption

# define number of pixels
config_list = np.sort(glob.glob(case_dir + "meas_setup_*.dat"))
nfiles = len(config_list)

# define iterable arguments for multiprocessing module
base_dir_iter = np.repeat(base_dir,nfiles)
libradtran_dir_iter = np.repeat(libradtran_dir,nfiles)
pix_id_iter = np.arange(1,nfiles+1,1)

# initialize multiprocessing module
with concurrent.futures.ProcessPoolExecutor(max_workers=n_core) as executor:
    
    # running retrieval l2m in parallel
    executor.map(l2m_wrapper, base_dir_iter, libradtran_dir_iter, pix_id_iter)

# print statement
print("Info         | ")
print("Info         |", dashed)
print("Info         | End of Level-2 Retrieval Module")
print("Info         |", dashed)