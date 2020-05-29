#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np

def read_uvspec_output(input_file, output_file):
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++
    # translate the output to desired radiative quantities
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++
    # open and read output file
    dummy = open(input_file,"rt")
    data = dummy.readlines()
    dummy.close()
    
    # define wavelength
    step = 3    # 1 umu = 3 steps, 2 umu = 4 steps
    idx_arr = np.arange(0,len(data),step)
    nlambda = len(idx_arr)

    # create array
    output_array = np.zeros(shape=(nlambda,3), dtype=float)

    # write output file
    file = open(output_file, "w")
    file.write(" lambda   Irr_dw            Rad_up \n")

    # define the format
    output_format = "%8.3f  %15.10e  %15.10e \n"
    
    # loop over lambda
    for i in range(nlambda):
        # select index
        idx = idx_arr[i]
        
        # parse output and pass them to the array (in Watts)
        output_array[i,0] = float(data[idx].split()[0])
        output_array[i,1] = float(data[idx].split()[1])/1E3+float(data[idx].split()[2])/1E3
        output_array[i,2] = float(data[idx+2].split()[2])/1E3
        
        # write output
        file.write(output_format %(output_array[i,0], output_array[i,1], output_array[i,2]))
    
    # close file    
    file.close()