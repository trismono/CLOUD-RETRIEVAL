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

# define input parameters
libradtran_dir = "/deos/trismonock/libRadtran-2.0.2/"
atm_profile = "afglms"
source = "solar"
solar_type = "kurudz_1.0nm"
rte_solver = "fdisort2"
mol_abs_parameter = "lowtran"
umu = "-1.0, 1.0"               # default setting = upward and nadir looking
wc_mode = "mie"                 # for liquid water cloud 
nstream = "16"                  # for better radiance accuracy, increase nstream (number of stream)
ncore = "3"                     # indeally, this should be number of cpu-1

# define filename
filename = input_dir + "fwd_config.dat"

# create input directory
try:
    os.makedirs(input_dir)
except OSError:
    pass

# create config file
config = ConfigParser()

# libradtran
config["LIBRADTRAN"] = {
    "libradtran_dir": "%s" %libradtran_dir,
    "atm_profile": "%s" %atm_profile,
    "source": "%s" %source,
    "solar_type": "%s" %solar_type,
    "rte_solver": "%s" %rte_solver,
    "mol_abs_parameter": "%s" %mol_abs_parameter,
    "wc_mode": "%s" %wc_mode,
    "umu": "%s" %umu,
    "nstream": "%s" %nstream,
    }

# number of core
config["CORE"] = {
    "n_core": "%s" %ncore
    }

# writing file
with open(filename, "w") as f:
    config.write(f)