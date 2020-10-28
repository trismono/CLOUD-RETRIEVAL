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
import time

# import local modules
from read_case import read_case
from uvspec_input import uvspec_input
from cloud_profile import cloud_profile
from read_uvspec_output import read_uvspec_output

def fwd_wrapper(base_dir, libradtran_dir, config_file, pix_id):
    # define process id
    pid = os.getpid()
    
    # define paths
    case_dir = os.path.realpath(base_dir + "../../01_CREATE_CASE/OUTPUT") + "/"
    input_dir = os.path.realpath(base_dir + "../INPUT") + "/"
    output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"
    
    # libtradtran internal database
    lib_data_dir = os.path.realpath(libradtran_dir + "data")   

    # print statement
    print("Info         | PID %i processing pixel number : %i" %(pid, pix_id))
    
    # define config file
    case_file = case_dir + "meas_setup_%i.dat" %pix_id

    # define time start
    time0 = time.perf_counter()
    
    # +++++++++++++++++++++
    # read and parse config
    # +++++++++++++++++++++
    # read config
    (sza,phi0,phi,zout,doy,albedo_file,ztop,zbase,reff,tau550,aerosol_season,\
     aerosol_haze,lambda0,lambda1) = read_case(case_file=case_file)
        
    # adjust format
    sza =  str(sza)
    phi0 = str(phi0)
    phi = str(phi)
    zout = str(zout)
    doy = str(doy)
    albedo_file = str(albedo_file)
    ztop = float(ztop)
    zbase = float(zbase)
    reff = float(reff)
    tau550 = float(tau550)
    aerosol_season = str(aerosol_season)
    aerosol_haze = str(aerosol_haze)
    lambda0 = str(lambda0)
    lambda1 = str(lambda1)

    # +++++++++++++++++++
    # build cloud profile
    # +++++++++++++++++++
    # define cloud file
    cloud_file = input_dir + "cloud_%.2f_%.2f.dat" %(tau550,reff)
    
    # define lwc (not important:: in the end, lwc will be modified as tau550 is set)
    lwc = 0.2 

    # create cloud profile    
    cloud_profile(cloud_file=cloud_file,ztop=ztop,zbase=zbase,lwc=lwc,reff=reff)

    # +++++++++++++++++++++    
    # generate uvspec input
    # +++++++++++++++++++++    
    # input filename
    input_file = input_dir + "meas_setup_%i.inp" %pix_id

    # uvspec input    
    uvspec_input(config_file=config_file,lib_data_dir=lib_data_dir,cloud_file=cloud_file,\
                 input_file=input_file,sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,\
                 albedo_file=albedo_file,lambda0=lambda0,lambda1=lambda1,tau550=tau550,\
                 aerosol_season=aerosol_season,aerosol_haze=aerosol_haze)

    # +++++++++++   
    # run unvspec
    # +++++++++++
    # define output file
    output_file = output_dir + "meas_setup_%i.out" %pix_id

    # verbose file 
    # verbose_file = output_dir + "meas_setup_%i_verbose.txt" %pix_id
    
    # define run script (use verbose only if you wan't it)
    # run_script = "(uvspec <" + input_file + "> " + output_file + ") >& " + verbose_file
    run_script = "uvspec <" + input_file + "> " + output_file

    # run uvspec    
    os.system(run_script)
    
    # ++++++++++++++++++++++    
    # reformat uvspec output 
    # ++++++++++++++++++++++
    # define output file
    output_file_fmt = output_dir + "meas_setup_%i.dat" %pix_id
    
    # reading and reformat uvspec output
    read_uvspec_output(input_file=output_file, output_file=output_file_fmt)
    
    # define end time
    time1 = time.perf_counter()
    
    # calculate execution time
    time_total = time1-time0
    
    # print statement
    print("Info         | PID %i elapsed time : %.1f sec" %(pid, time_total))
    
    