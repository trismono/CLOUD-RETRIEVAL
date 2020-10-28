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

# limit openblas thread number
os.environ["OPENBLAS_NUM_THREADS"] = "3"

import numpy as np
import glob
from configparser import ConfigParser

# import local routines
from l1b_convolution import l1b_convolution
from l1b_noise import l1b_noise
from l1b_interface import l1b_interface

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define paths
data_dir = os.path.realpath(base_dir + "../../DATA") + "/"
fwd_output_dir = os.path.realpath(base_dir + "../../02_SIM_FWD/OUTPUT") + "/"
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"
output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"

# define config file
config_file = input_dir + "l1b_config.dat"

# create dir
try:
    os.makedirs(output_dir)
except OSError:
    pass

# open config file
config = ConfigParser()
config.read(config_file)

# parse config
snr_wavelength = config['SNR']['wavelength']
snr_smart = config['SNR']['snr_smart']
fwhm = config['ISRF']['fwhm']
nisrf_grid = config['ISRF']['nisrf_grid']
lambda_offset = config['ISRF']['lambda_offset']

# print statement
dashed = 60 * "+"
print("Info         |", dashed)
print("Info         | End-to-End Simulator Level-1B Module")
print("Info         | Author       : Trismono C. Krisna")
print("Info         | Institution  : SRON Netherlands Institute for Space Research")
print("Info         | Version      : v2.2")
print("Info         | Date         : 27 October 2020")
print("Info         |", dashed)
print("Info         | ")
print("Info         | Wavelengths of SNR :", snr_wavelength)
print("Info         | SNR :", snr_smart)
print("Info         | FWHM of ISRF :", fwhm)
print("Info         | Number of ISRF grid :", nisrf_grid)
print("Info         | Offset from ISRF center :", lambda_offset)
print("Info         | ")

# +++++++++++++++++++++++
# read spectral grid data
# +++++++++++++++++++++++
# open and read data
grid_vnir = np.loadtxt(data_dir + "wl_pixel_Iup_vnir.dat", skiprows=2)
grid_swir = np.loadtxt(data_dir + "wl_pixel_Iup_swir.dat", skiprows=2)

# define wavelengths
lambda_vnir = grid_vnir[:,1]
lambda_swir = np.flipud(grid_swir[:,1])

# define dimension of spectra region: vnir and swir
nspectra_region = 2

# +++++++++++++++++++++++++++++++
# running l1b module: convolution
# +++++++++++++++++++++++++++++++
# search measurement setup files
file_list = np.sort(glob.glob(fwd_output_dir + "meas_setup_*.dat"))

# define dimension
nfiles = len(file_list) 

# loop over files
for i in range(nfiles):
    # print statement
    print("Info         | Executing convolution module for pixel id:", i+1)
    
    # define input file
    input_file = fwd_output_dir + "meas_setup_%i.dat" %(i+1)

    # loop over spectra region  
    for j in range(nspectra_region):
        # for vnir
        if j == 0:
            # define output file
            output_file = output_dir + "meas_setup_%i_vnir.dat" %(i+1)
            
            # run l1b module
            l1b_convolution(config_file=config_file, input_file=input_file,spectral_grid=lambda_vnir, output_file=output_file)

        # for swir
        elif j ==1:
            # define output file
            output_file = output_dir + "meas_setup_%i_swir.dat" %(i+1)
            
            # run l1b module
            l1b_convolution(config_file=config_file, input_file=input_file,spectral_grid=lambda_swir, output_file=output_file)            

# ++++++++++++++++++++++++++++++++++++++++++++++++
# running l1b module: radiance noise and interface
# ++++++++++++++++++++++++++++++++++++++++++++++++
# note: noise realization is based on MODIS SNR
# assuming that SMART noise is two timer larger than MODIS
# therefore SMART SNR is half of MODIS SNR

# loop over files
for i in range(nfiles):
    # print statement
    print("Info         | Executing noise module for pixel id:", i+1)
    
    # loop over spectra region
    for j in range(nspectra_region):
        # for vnir
        if j == 0:
            # define filename
            input_file = output_dir + "meas_setup_%i_vnir.dat" %(i+1)
            
            # define l1b filename
            output_file = output_dir + "radiance_meas_%i_vnir.dat" %(i+1)   
            
        # for swir            
        elif j == 1:
            # define filename
            input_file = output_dir + "meas_setup_%i_swir.dat" %(i+1)

            # define l1b filename
            output_file = output_dir + "radiance_meas_%i_swir.dat" %(i+1)
            
        # open and read spectra
        spectra_data = np.loadtxt(input_file, skiprows=1)

        # define radiance and wavelengths
        wavelength = spectra_data[:,0]
        rad_up = spectra_data[:,4]
        
        # run l1b noise
        (noisy_spectra, noise) = l1b_noise(config_file=config_file, wavelength=wavelength, spectra=rad_up)
        
        # write l1b interface
        l1b_interface(wavelength=wavelength, noisy_spectra=noisy_spectra, pure_spectra=rad_up, noise=noise, output_file=output_file)

# print statement
print("Info         | ")
print("Info         |", dashed)
print("Info         | End of Simulator Level-1B Module")
print("Info         |", dashed)