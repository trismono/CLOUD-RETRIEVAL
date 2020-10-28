#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 08:34:01 2020

@author: trismonok
"""

def write_spectra_to_ascii(filename, data):
    # define dimension
    nspectra = len(data)
    
    # open file
    file = open(filename, "w")
    
    # write header
    file.write("%8s %18s %18s\n" %("lambda", "Irr_dw", "Rad_up"))
    
    # define the format
    output_format = "%8.3f %18.10e %18.10e\n"
    
    # loop over spectra
    for i in range(nspectra):
        # write output
        file.write(output_format %(data[i,0], data[i,1], data[i,2]))
    
    # close file
    file.close()