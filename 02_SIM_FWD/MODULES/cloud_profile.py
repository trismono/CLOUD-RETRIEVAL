#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

def cloud_profile(cloud_file,ztop,zbase,lwc,reff):    
    # open file
    file = open(cloud_file, "w")
    
    # writing parameter    
    file.write("#  z(km)  LWC(g/m3)  r_eff(micron)\n")
    file.write("%6.3f     %6.3f         %6.3f\n"    %(ztop, 0, 0))           
    file.write("%6.3f     %6.3f         %6.3f"      %(zbase, lwc, reff))
    
    # close file
    file.close()   