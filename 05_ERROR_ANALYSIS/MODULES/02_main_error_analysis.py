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
from numpy.linalg import inv
from configparser import ConfigParser

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define dirs
kmat_dir = os.path.realpath(base_dir + "../../04_SIM_RET/OUTPUT") + "/"
input_dir = os.path.realpath(base_dir + "../INPUT") + "/"
output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"

# create dir
try:
    os.makedirs(output_dir)
except OSError:
    pass

# print statement
dashed = 60 * "+"
print("Info         |", dashed)
print("Info         | End-to-End Simulator Error Analysis Module")
print("Info         | Author       : Trismono C. Krisna")
print("Info         | Institution  : SRON Netherlands Institute for Space Research")
print("Info         | Version      : v2.2")
print("Info         | Date         : 27 October 2020")
print("Info         |", dashed)

# ++++++++++++
# read setting
# ++++++++++++
# define config file
config_file = input_dir + "error_config.dat"

# open config file
config = ConfigParser()
config.read(config_file)

# parse config
weighting_tau = float(config['WEIGHTING']['tau'])
weighting_reff = float(config['WEIGHTING']['reff'])

# print statement
print("Info         | ")
print("Info         | Variance of tau :", weighting_tau)
print("Info         | Variance of reff :", weighting_reff, "micron")
print("Info         | ") 

# define cases
stat_list = np.sort(glob.glob(kmat_dir + "statistics*.dat"))
nfiles = len(stat_list)

# loop over pixels
for i in range(nfiles):
    # define init string
    kmat_list = glob.glob(kmat_dir + "kmat_" + str(i+1).zfill(5) + "_*.dat")
    
    # define number of iteration
    nkmat = len(kmat_list)
    
    # define filename
    kmat_file = kmat_dir + "kmat_" + str(i+1).zfill(5) + "_" + str(nkmat).zfill(2) + ".dat"
    
    # open and read data
    data = np.loadtxt(kmat_file, skiprows=1)
    
    wavelength = data[:,0]
    y_meas = data[:,1]
    meas_error = data[:,2]
    y_mod = data[:,3]
    kmat = data[:,4:]
    
    # define measurement covariance sy
    sy = np.zeros(shape=(len(wavelength),len(wavelength)), dtype=float)
    np.fill_diagonal(sy,meas_error)
        
    # define prior covariance sa :: should consistent with the retrieval setting
    weighting = np.array([weighting_tau,weighting_reff], dtype=float)
    sa = np.zeros(shape=(len(weighting),len(weighting)), dtype=float)
    np.fill_diagonal(sa,weighting**2)
    
    # ++++++++++++++
    # error analysis
    # ++++++++++++++
    # Covariance matrix of the solution
    # Sx = (Kt Sy-1 K + Sa-1)-1
    sx = inv(np.matmul(kmat.T,np.matmul(inv(sy),kmat)) + inv(sa))
    
    # Gain-matrix
    # G = (Kt Sy-1 K + Sa-1)-1 Kt Sy-1 = Sx Kt Sy-1
    G = np.matmul(sx,np.matmul(kmat.T,inv(sy)))
    
    # Averaging kernel
    # A = G K
    A = np.matmul(G,kmat)  
    
    # Smooting error
    # Sx_s = (I - A) Sa (I - A)t
    sx_smooth = np.matmul(np.identity(len(sa))-A,np.matmul(sa,(np.identity(len(sa))-A).T))  
    
    # Noise error
    # Sx_n = G Se Gt
    sx_noise = np.matmul(G,np.matmul(sy,np.transpose(G)))
    
    # +++++++++++++
    # final results
    # +++++++++++++
    # parse
    prior_error = weighting
    total_error = np.sqrt(np.diag(sx))
    noise_error = np.sqrt(np.diag(sx_noise))
    smoothing_error = np.sqrt(np.diag(sx_smooth))
    
    # print results
    print("Info         | Pixel ID %i :: Prior error tau = %.2f and reff = %.2f micron" %(i+1, weighting[0], weighting[1]))
    print("Info         | Pixel ID %i :: Noise error tau = %.2f and reff = %.2f micron" %(i+1, noise_error[0], noise_error[1]))
    print("Info         | Pixel ID %i :: Smoothing error tau = %.2f and reff = %.2f micron" %(i+1, smoothing_error[0], smoothing_error[1]))
    print("Info         | Pixel ID %i :: Total error tau = %.2f and reff = %.2f micron" %(i+1, total_error[0], total_error[1]))
    
    # ++++++++++++++
    # writing output
    # ++++++++++++++
    # join data
    data_dummy = np.array([prior_error, noise_error, smoothing_error, total_error]).T
    
    # define filename
    out_file = output_dir + "error_statistics_" + str(i+1).zfill(5) + ".dat"

    # write outputs    
    with open(out_file, 'w+') as fid:
        # make header
        fid.writelines("+++++++++++++++++++++\n")
        fid.writelines("Linear Error Analysis\n")
        fid.writelines("+++++++++++++++++++++\n")
        fid.writelines("Col 1 = Prior_Error\n")
        fid.writelines("Col 2 = Noise_Error\n")
        fid.writelines("Col 3 = Smoothing_Error\n")
        fid.writelines("Col 4 = Total_Error\n")
        fid.writelines("Row 1 = Optical thickness\n")
        fid.writelines("Row 2 = Effective radius [micron]\n")

        # write data
        np.savetxt(fid, data_dummy, fmt=['%17.10e','%17.10e','%17.10e','%17.10e'])
    
# print statement
print("Info         | ")
print("Info         |", dashed)
print("Info         | End of Error Analysis Module")
print("Info         |", dashed)