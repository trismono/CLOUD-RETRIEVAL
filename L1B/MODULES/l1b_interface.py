#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

def l1b_interface(wavelength, noisy_spectra, pure_spectra, noise, output_file):
    # define dimension
    nspectra = len(wavelength)
    
    # write output file
    file = open(output_file, "w")
    file.write(" lambda     Rad_up_noisy   Rad_up_pure    Rad_up_noise\n")

    # define the format
    output_format = "%8.3f  %13.5e  %13.5e  %13.5e\n"
    
    # loop over spectra
    for i in range(nspectra):
        # write output
        file.write(output_format %(wavelength[i], noisy_spectra[i], \
                                   pure_spectra[i], noise[i]))

    # close file
    file.close()
    
    
    