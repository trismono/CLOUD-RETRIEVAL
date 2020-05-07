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

# import local routines
from l1b_convolution import l1b_convolution
from l1b_noise import l1b_noise
from l1b_interface import l1b_interface

# get and define the working directory
dir_path = os.getcwd() + "/"
os.chdir(dir_path)

# define paths
data_dir = os.path.realpath(dir_path + "../../DATA") + "/"
input_dir = os.path.realpath(dir_path + "../../SIM_FWD/INOUT") + "/"
output_dir = os.path.realpath(dir_path + "../OUTPUT") + "/"

# create dir
try:
    os.makedirs(output_dir)
except OSError:
    pass

# +++++++++++++++++++++++
# read spectral grid data
# +++++++++++++++++++++++
# open and read data
grid_vnir = np.loadtxt(data_dir + "wl_pixel_Iup_vnir.dat", skiprows=2)
grid_swir = np.loadtxt(data_dir + "wl_pixel_Iup_swir.dat", skiprows=2)

# define wavelengths
lambda_vnir = grid_vnir[:,1]
lambda_swir = np.flipud(grid_swir[:,1])

# define dimension spectral grid
nlambda_grid = 2

# +++++++++++++++++++++++++++++++
# running l1b module: convolution
# +++++++++++++++++++++++++++++++
# search measurement setup files
file_list = glob.glob(input_dir + "meas_setup_*.dat")
nfiles = len(file_list) 

# loop over files
for i in range(nfiles):
    # define input file
    input_file = input_dir + "meas_setup_%i.dat" %(i+1)
    
    for j in range(nlambda_grid):
        if j == 0:
            # define output file
            output_file = output_dir + "meas_setup_%i_vnir.dat" %(i+1)
            
            # run l1b module
            l1b_convolution(input_file=input_file,spectral_grid=lambda_vnir,\
                            output_file=output_file)
        elif j ==1:
            # define output file
            output_file = output_dir + "meas_setup_%i_swir.dat" %(i+1)
            
            # run l1b module
            l1b_convolution(input_file=input_file,spectral_grid=lambda_swir,\
                            output_file=output_file)            

# ++++++++++++++++++++++++++++++++++++++++++++++++
# running l1b module: radiance noise and interface
# ++++++++++++++++++++++++++++++++++++++++++++++++
"""
Note from Author:
- Noise realization is based on MODIS SNR (simplified model)
- Assuming that SMART noise is two timer larger than MODIS
- Therefore SMART SNR is half of MODIS SNR
- Of course you can define your own SNR and noise
"""

# loop over files
for i in range(nfiles):
    # loop over lambda grid
    for j in range(nlambda_grid):
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
        (noisy_spectra, noise) = l1b_noise(wavelength=wavelength, spectra=rad_up)
        
        # write l1b interface
        l1b_interface(wavelength=wavelength, noisy_spectra=noisy_spectra, \
                      pure_spectra=rad_up, noise=noise, output_file=output_file)