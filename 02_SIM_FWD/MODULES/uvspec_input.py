#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

from configparser import ConfigParser

def uvspec_input(config_file,lib_data_dir,cloud_file,input_file,sza,phi0,phi,zout,doy,\
                 albedo_file,lambda0,lambda1,tau550,aerosol_season,aerosol_haze):
    
    # ++++++++++++++++
    # read config file
    # ++++++++++++++++    
    # open config file
    config = ConfigParser()
    config.read(config_file)
    
    # parse config
    atm_profile = config['LIBRADTRAN']['atm_profile']
    source = config['LIBRADTRAN']['source']
    solar_type = config['LIBRADTRAN']['solar_type']
    rte_solver = config['LIBRADTRAN']['rte_solver']
    mol_abs_parameter = config['LIBRADTRAN']['mol_abs_parameter']
    wc_mode = config['LIBRADTRAN']['wc_mode']
    umu = (config['LIBRADTRAN']['umu']).split(",")
    nstream = config['LIBRADTRAN']['nstream']
    
    # +++++++++++++++++
    # create input file
    # +++++++++++++++++
    # open file
    file = open(input_file, "w") 

    # write input parameters       
    file.write("data_files_path " + lib_data_dir +"/\n")
    file.write("atmosphere_file " + lib_data_dir + "/atmmod/" + atm_profile + ".dat\n")
    file.write("source " + source + " " + lib_data_dir + "/solar_flux/" + solar_type + ".dat\n")
    file.write("albedo_file " + albedo_file + "\n")
    file.write("day_of_year " + doy + "\n")
    file.write("sza " + sza + "\n")
    file.write("phi0 " + phi0 + "\n")    
    file.write("phi " + phi + "\n")
    file.write("umu " + umu[0] + umu[1] + " \n")
    file.write("zout " + zout + "\n")
    file.write("rte_solver " + rte_solver + "\n")
    file.write("mol_abs_param " + mol_abs_parameter + "\n")
    file.write("number_of_streams " + nstream + "\n")
    file.write("wavelength "  + lambda0 + " " + lambda1 + "\n")
    file.write("aerosol_default\n")
    file.write("aerosol_haze " + aerosol_haze + "\n")
    file.write("aerosol_season " + aerosol_season + "\n")
    file.write("wc_file 1D " + cloud_file + "\n") 
    file.write("wc_properties " + wc_mode + " interpolate\n")
    file.write("wc_modify tau550 set %.2f\n" %tau550) 	     
    file.write("quiet")         # quiet
    # file.write("verbose")       # print verbose file (useful for post analysis)

    # close file
    file.close()