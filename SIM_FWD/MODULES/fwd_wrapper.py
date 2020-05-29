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
from read_config import read_config
from uvspec_input import uvspec_input
from cloud_profile import cloud_profile
from read_uvspec_output import read_uvspec_output

def fwd_wrapper(base_dir, pix_id):
    # define paths
    case_dir = os.path.realpath(base_dir + "../../CREATE_CASE/OUTPUT") + "/"
    inout_dir = os.path.realpath(base_dir + "../INOUT") + "/"
    cloud_dir = os.path.realpath(base_dir + "../CLOUD") + "/"
    lib_dir = "/deos/trismonock/libRadtran-2.0.2/"      # !! this must be adjusted !!
    lib_data_dir = os.path.realpath(lib_dir + "data")   # libtradtran internal database

    # print statement
    print("Info     | Processing pixel number :", pix_id)
    
    # define config file
    config_file = case_dir + "meas_setup_%i.dat" %pix_id
    
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
    
    # define lwc (not important:: in the end, lwc will be modified as tau550 is set)
    lwc = 0.2 

    # create cloud profile    
    cloud_profile(cloud_file=cloud_file,ztop=float(ztop),zbase=float(zbase),\
                  lwc=lwc,reff=float(reff))

    # +++++++++++++++++++++    
    # generate uvspec input
    # +++++++++++++++++++++    
    # input filename
    input_file = inout_dir + "meas_setup_%i.inp" %pix_id

    # uvspec input    
    uvspec_input(lib_data_dir=lib_data_dir,cloud_file=cloud_file,input_file=input_file,\
                 sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,albedo_file=albedo_file,\
                 lambda0=lambda0,lambda1=lambda1,tau550=tau550,aerosol_season=aerosol_season,\
                 aerosol_haze=aerosol_haze)
    
    # +++++++++++   
    # run unvspec
    # +++++++++++
    # define output file
    output_file = inout_dir + "meas_setup_%i.out" %pix_id

    # verbose file 
    verbose_file = inout_dir + "meas_setup_%i_verbose.txt" %pix_id
    
    # define run script (use verbose only if you wan't it)
    # run_script = "(uvspec <" + input_file + "> " + output_file + ") >& " + verbose_file
    run_script = "uvspec <" + input_file + "> " + output_file

    # run uvspec    
    os.system(run_script)
    
    # ++++++++++++++++++++++    
    # reformat uvspec output 
    # ++++++++++++++++++++++
    # define output file
    output_file_fmt = inout_dir + "meas_setup_%i.dat" %pix_id
    
    # reading and reformat uvspec output
    read_uvspec_output(input_file=output_file, output_file=output_file_fmt)
    
    
    