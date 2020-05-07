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
import glob

# import local routines
from read_config import read_config
from uvspec_input import uvspec_input
from cloud_profile import cloud_profile
from read_uvspec_output import read_uvspec_output

# get and define the working directory
dir_path = os.getcwd() + "/"
os.chdir(dir_path)

# define paths
case_dir = os.path.realpath(dir_path + "../../CREATE_CASE/OUTPUT") + "/"
data_dir = os.path.realpath(dir_path + "../../DATA") + "/"
inout_dir = os.path.realpath(dir_path + "../INOUT") + "/"
cloud_dir = os.path.realpath(dir_path + "../CLOUD") + "/"
lib_dir = "/deos/trismonock/libRadtran-2.0.2/"      # !! this must be adjusted !!
lib_data_dir = os.path.realpath(lib_dir + "data")   # libtradtran internal database

# create input-output directory
try:
    os.makedirs(inout_dir)
    os.makedirs(cloud_dir)
except OSError:
    pass

# ++++++++++++++++++++++++++++
# define cases to be simulated
# ++++++++++++++++++++++++++++
# search measurement setup files
config_list = glob.glob(case_dir + "meas_setup_*.dat")
nfiles = len(config_list) 

# loop over pixels
for i in range(nfiles):
    # print statement
    print("Info     | Processing pixel number :", i+1)
    
    # define config file
    config_file = case_dir + "meas_setup_%i.dat" %(i+1)
    
    # +++++++++++++++++++++
    # read and parse config
    # +++++++++++++++++++++
    (sza,phi0,phi,zout,doy,albedo_file,ztop,zbase,reff,tau550,aerosol_season,\
     aerosol_haze,lambda0,lambda1) = read_config(config_file=config_file)

    # +++++++++++++++++++
    # build cloud profile
    # +++++++++++++++++++
    # define cloud file
    cloud_file = cloud_dir + "cloud_%.2f_%.2f.dat" %(float(tau550),float(reff))
    
    # define lwc
    lwc = 0.2   # in the end, lwc will be modified as tau550 is set

    # create cloud profile    
    cloud_profile(cloud_file=cloud_file,ztop=float(ztop),zbase=float(zbase),\
                  lwc=lwc,reff=float(reff))

    # +++++++++++++++++++++    
    # generate uvspec input
    # +++++++++++++++++++++    
    # input filename
    input_file = inout_dir + "meas_setup_%i.inp" %(i+1)

    # uvspec input    
    uvspec_input(lib_data_dir=lib_data_dir,cloud_file=cloud_file,input_file=input_file,\
                 sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,albedo_file=albedo_file,\
                 lambda0=lambda0,lambda1=lambda1,tau550=tau550,aerosol_season=aerosol_season,\
                 aerosol_haze=aerosol_haze)
    
    # +++++++++++   
    # run unvspec
    # +++++++++++
    # define output file
    output_file = inout_dir + "meas_setup_%i.out" %(i+1)

    # verbose file 
    verbose_file = inout_dir + "meas_setup_%i_verbose.txt" %(i+1)
    
    # define run script (use verbose only if you wan't it)
#    run_script = "(uvspec <" + input_file + "> " + output_file + ") >& " + verbose_file
    run_script = "uvspec <" + input_file + "> " + output_file

    # run uvspec    
    os.system(run_script)
    
    # ++++++++++++++++++++++    
    # reformat uvspec output 
    # ++++++++++++++++++++++
    # define output file
    output_file_fmt = inout_dir + "meas_setup_%i.dat" %(i+1)
    
    # reading and reformat uvspec output
    read_uvspec_output(input_file=output_file, output_file=output_file_fmt)
    
    
    