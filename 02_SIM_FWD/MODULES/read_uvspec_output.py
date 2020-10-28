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
    step = 4
    idx_arr = np.arange(0,len(data),step)
    nlambda = len(idx_arr)

    # create array
    output_array = np.zeros(shape=(nlambda,11), dtype=float)

    # write output file
    file = open(output_file, "w")
    file.write(" lambda   Irr_dw       Irr_up       Rad_dw       Rad_up       Act_dw       Act_up       Irr_dw_dir   Irr_up_diff  Act_dw_dir   Act_up_diff\n")

    # define the format
    output_format = "%8.3f  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n"
    
    # loop over lambda
    for i in range(nlambda):
        # select index
        idx = idx_arr[i]
        
        # parse output and pass them to the array (in Watts)
        output_array[i,0] = float(data[idx].split()[0])
        output_array[i,1] = float(data[idx].split()[1])/1E3+float(data[idx].split()[2])/1E3
        output_array[i,2] = float(data[idx].split()[3])/1E3
        output_array[i,3] = float(data[idx+2].split()[2])/1E3
        output_array[i,4] = float(data[idx+3].split()[2])/1E3
        output_array[i,5] = float(data[idx].split()[4])/1E3+float(data[idx].split()[5])/1E3             
        output_array[i,6] = float(data[idx].split()[6])/1E3
        output_array[i,7] = float(data[idx].split()[1])/1E3
        output_array[i,8] = float(data[idx].split()[2])/1E3
        output_array[i,9] = float(data[idx].split()[4])/1E3
        output_array[i,10] = float(data[idx].split()[5])/1E3
        
        # write output
        file.write(output_format %(output_array[i,0], output_array[i,1], output_array[i,2], \
                                   output_array[i,3], output_array[i,4], output_array[i,5], \
                                   output_array[i,6], output_array[i,7], output_array[i,8], \
                                   output_array[i,9], output_array[i,10]))
    # close file
    file.close()