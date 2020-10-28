#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import numpy as np

def read_uvspec_output(input_file):
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
    
    # loop over lambda
    for i in range(nlambda):
        # select index
        idx = idx_arr[i]
        
        # parse output and pass them to the array (in Watts)
        output_array[i,0] = float(data[idx].split()[0])
        output_array[i,1] = float(data[idx].split()[1])/1E3+float(data[idx].split()[2])/1E3
        output_array[i,2] = float(data[idx+2].split()[2])/1E3
        
    # return output
    return(output_array)
