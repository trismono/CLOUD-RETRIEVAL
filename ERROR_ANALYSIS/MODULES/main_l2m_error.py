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
from numpy.linalg import inv

# get and define the working directory
base_dir = os.getcwd() + "/"
os.chdir(base_dir)

# define dirs
input_dir = os.path.realpath(base_dir + "../../SIM_RET/OUTPUT") + "/"
output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"

# create dir
try:
    os.makedirs(output_dir)
except OSError:
    pass

# define cases
stat_list = glob.glob(input_dir + "statistics*.dat")
stat_list = np.sort(stat_list)
nfiles = len(stat_list)

# loop over pixels
for i in range(nfiles):
    # print statement
    print("Info    | Processing pixel id: ", i+1)

    # define init string
    kmat_list = glob.glob(input_dir + "kmat_" + str(i+1).zfill(5) + "_*.dat")
    
    # define number of iteration
    nkmat = len(kmat_list)
    
    # define filename
    kmat_file = input_dir + "kmat_" + str(i+1).zfill(5) + "_" + str(nkmat).zfill(2) + ".dat"
    
    # open and read data
    data = np.loadtxt(kmat_file)
    
    wavelength = data[:,0]
    y_meas = data[:,1]
    meas_error = data[:,2]
    y_mod = data[:,3]
    kmat = data[:,4:]
    
    # define measurement covariance sy
    sy = np.zeros(shape=(len(wavelength),len(wavelength)), dtype=float)
    np.fill_diagonal(sy,meas_error)
        
    # define prior covariance sa :: should consistent with the retrieval setting
    weighting = np.array([16.0,20.0], dtype=float)
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
    prior_error = weighting
    total_error = np.sqrt(np.diag(sx))
    noise_error = np.sqrt(np.diag(sx_noise))
    smoothing_error = np.sqrt(np.diag(sx_smooth))
    
    # ++++++++++++++
    # writing output
    # ++++++++++++++
    # join data
    data_dummy = np.array([prior_error, noise_error, smoothing_error, total_error]).T
    
    # define filename
    out_file = output_dir + "ret_error_" + str(i+1).zfill(5) + ".dat"

    # write outputs    
    with open(out_file, 'w+') as fid:
        # write a header variable names and units, and data
        fid.writelines(["%17s %17s %17s %17s \n" %("Prior error","Noise error","Smoothing error","Total error")])
        np.savetxt(fid, data_dummy, fmt=['%17.10e','%17.10e','%17.10e','%17.10e'])
    
