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
    
    # open file
    file = open(output_file, "w")
    
    # write header
    file.write("%8s %15s %15s %15s\n" %("lambda", "Rad_up_noisy", "Rad_up_pure", "Rad_up_noise"))

    # define format
    output_format = "%8.3f %15.5e %15.5e %15.5e\n"
    
    # loop over spectra
    for i in range(nspectra):
        # write output
        file.write(output_format %(wavelength[i], noisy_spectra[i], pure_spectra[i], noise[i]))

    # close file
    file.close()
    
    
    