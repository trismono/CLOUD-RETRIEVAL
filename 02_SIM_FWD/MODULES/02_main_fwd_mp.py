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

# import local routines
from fwd_wrapper import fwd_wrapper

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define paths
case_dir = os.path.realpath(base_dir + "../../01_CREATE_CASE/OUTPUT") + "/"
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"
output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"

# create input-output directory
try:
    os.makedirs(input_dir)
    os.makedirs(output_dir)
except OSError:
    pass

# print statement
dashed = 60 * "+"
print("Info         |", dashed)
print("Info         | End-to-End Simulator Forward Simulation")
print("Info         | Author       : Trismono C. Krisna")
print("Info         | Institution  : SRON Netherlands Institute for Space Research")
print("Info         | Version      : v2.2")
print("Info         | Date         : 27 October 2020")
print("Info         |", dashed)

# ++++++++++++
# read setting
# ++++++++++++
# define config file
config_file = input_dir + "fwd_config.dat"

# open config file
config = ConfigParser()
config.read(config_file)

# parse config
libradtran_dir = config['LIBRADTRAN']['LIBRADTRAN_DIR']
atm_profile = config['LIBRADTRAN']['atm_profile']
source = config['LIBRADTRAN']['source']
solar_type = config['LIBRADTRAN']['solar_type']
rte_solver = config['LIBRADTRAN']['rte_solver']
mol_abs_parameter = config['LIBRADTRAN']['mol_abs_parameter']
wc_mode = config['LIBRADTRAN']['wc_mode']
umu = config['LIBRADTRAN']['umu']
nstream = config['LIBRADTRAN']['nstream']
n_core = int(config['CORE']['n_core'])

# print statement
print("Info         | ")
print("Info         | LibRadtran dir :", libradtran_dir)
print("Info         | Atmospheric profile :", atm_profile)
print("Info         | Radiation source :", source)
print("Info         | Solar type :", solar_type)
print("Info         | Radiative transfer solver :", rte_solver)
print("Info         | Molecular absorption :", mol_abs_parameter)
print("Info         | Scattering properties :", wc_mode)
print("Info         | Cosinus of viewing angle :", umu)
print("Info         | Number of core :", n_core)
print("Info         | ")

# ++++++++++++++++++++++
# multiprocessing module
# ++++++++++++++++++++++
# define number of process :: change this number accordingly
# when running the code in the local machine, optimum performance
# is normally n_proc = number of CPU (core) - 1
# when running the code on the server, please beware of
# processes by others and memory consumption

# define number of pixels
config_list = np.sort(glob.glob(case_dir + "meas_setup_*.dat"))
nfiles = len(config_list)

# define iterable arguments for multiprocessing module
base_dir_it = np.repeat(base_dir,nfiles)
libradtran_dir_it = np.repeat(libradtran_dir,nfiles)
config_file_it = np.repeat(config_file,nfiles)
pix_id_it = np.arange(1,nfiles+1,1)

# initialize multiprocessing module
with concurrent.futures.ProcessPoolExecutor(max_workers=n_core) as executor:
    
    # running retrieval l2m in parallel
    executor.map(fwd_wrapper, base_dir_it, libradtran_dir_it, config_file_it, pix_id_it) 

# print statement
print("Info         | ")
print("Info         |", dashed)
print("Info         | End of Simulator Forward Simulation")
print("Info         |", dashed)