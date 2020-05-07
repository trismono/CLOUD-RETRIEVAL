#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

import os, numpy as np, random
from read_config import read_ret_config, read_geometry_config
from uvspec_input import uvspec_input
from cloud_profile import cloud_profile
from convolution import l1b_convolution
from read_uvspec_output import read_uvspec_output
from kmat_module import kmat_module
from l2m_retrieval import l2m_retrieval
from diagnostic_module import diagnostic
from numpy.linalg import inv

def l2m_wrapper(base_dir,pix_id):     
    # ++++++++++++
    # define paths
    # ++++++++++++
    inp_dir = os.path.realpath(base_dir + "../INPUT") + "/"
    dummy_dir = os.path.realpath(base_dir + "../DUMMY") + "/"       # the size will be huge
    out_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"        # the size will be huge
    case_dir = os.path.realpath(base_dir + "../../CREATE_CASE/OUTPUT") + "/"
    l1b_dir = os.path.realpath(base_dir + "../../L1B/OUTPUT") + "/"
    lib_dir = "/deos/trismonock/libRadtran-2.0.2/"                  # !! this must be adjusted locally !!
    lib_data_dir = os.path.realpath(lib_dir + "data")               # libtradtran internal database
    
    # +++++++++++++++++++++
    # read retrieval config
    # +++++++++++++++++++++
    # define config file
    ret_config_file = inp_dir + "retrieval_setting.ini"
    
    # read and parse config
    (albedo_file,ztop,zbase,tau0,reff0,weighting_tau0,weighting_reff0,aerosol_season,\
     aerosol_haze,lambda0_fwd,lambda1_fwd,lambda0_ret_vnir,lambda1_ret_vnir,lambda0_ret_swir,\
     lambda1_ret_swir,flag_noise,rad_noise,iter_max) = read_ret_config(config_file=ret_config_file)
    
    # +++++++++++++++++++++++++
    # read measurement geometry
    # +++++++++++++++++++++++++
    # define geometry filename
    geometry_config_file = case_dir + "meas_setup_%i.dat" %pix_id
    
    # read measurement geometry
    (sza,phi0,phi,zout,doy) = read_geometry_config(config_file=geometry_config_file)
    
    # ++++++++++++++++
    # read measurement        
    # ++++++++++++++++
    # define filenames
    meas_vnir_file = l1b_dir + "radiance_meas_%i_vnir.dat" %pix_id
    meas_swir_file = l1b_dir + "radiance_meas_%i_swir.dat" %pix_id

    # read measurements        
    meas_vnir_data = np.loadtxt(meas_vnir_file, skiprows=1)
    meas_swir_data = np.loadtxt(meas_swir_file, skiprows=1)
    
    # read data
    lambda_meas_vnir = meas_vnir_data[:,0]
    lambda_meas_swir = meas_swir_data[:,0]
    
    # with noise or without noise
    if flag_noise == "1":
        rad_meas_vnir = meas_vnir_data[:,1]
        rad_meas_swir = meas_swir_data[:,1]
    elif flag_noise == "0":
        rad_meas_vnir = meas_vnir_data[:,2]
        rad_meas_swir = meas_swir_data[:,2]  
    else:
        print("Info     | Noise flag is not recognized")

    # ++++++++++++++++++++    
    # prepare measurements
    # ++++++++++++++++++++
    # find index
    idx_vnir = np.argwhere(np.logical_and(lambda_meas_vnir >= float(lambda0_ret_vnir), lambda_meas_vnir <= float(lambda1_ret_vnir)))
    idx_swir = np.argwhere(np.logical_and(lambda_meas_swir >= float(lambda0_ret_swir), lambda_meas_swir <= float(lambda1_ret_swir)))
    
    # arange measurement
    lambda_meas = np.append(lambda_meas_vnir[idx_vnir[:,0]],lambda_meas_swir[idx_swir[:,0]])
    rad_meas = np.append(rad_meas_vnir[idx_vnir[:,0]],rad_meas_swir[idx_swir[:,0]])
    
    # +++++++++++++++++++++++++++++++++
    # prepare covariance and init guess
    # +++++++++++++++++++++++++++++++++
    # measurement
    nmeas = len(rad_meas)
    sy = np.zeros(shape=(nmeas,nmeas), dtype=float)
    np.fill_diagonal(sy,(float(rad_noise)*rad_meas)**2)
    
    # state vector
    nstate = 2  # tau and reff
    sa = np.zeros(shape=(nstate,nstate), dtype=float)
    np.fill_diagonal(sa,np.array([float(weighting_tau0)**2,float(weighting_reff0)**2]))
    
    # define a priori
    xa = np.array([float(tau0),float(reff0)])
    
    # define initial guess x0 = xa
    x0 = xa     

    # +++++++++    
    # retrieval
    # +++++++++
    # define array
    xi_arr = np.zeros(shape=(int(iter_max)+1,2))
    chi2_arr = np.zeros(int(iter_max)+1)
    
    # assign x0 to xi
    xi_arr[0,:] = x0
    
    # ++++++++++++++++++
    # boundary condition
    # ++++++++++++++++++
    tau_min = 1
    tau_max = 50
    reff_min = 4        # lower boundary of mie scattering properties
    reff_max = 24       # upper boundary of mie scattering properties
    
    # start iteration : the number of iteration is set in the retrieval setting
    i = 0
    while i < int(iter_max):
        # print statement
        print("Info     | Pixel ID %i :: Iteration %i " %(pix_id,i+1))

        # +++++++++++++++++++++        
        # running forward model
        # +++++++++++++++++++++
        # define perturbation coefficient in percent
        c_ptau = np.array([0.06,0.04,0.02])
        c_ptau = np.array([0.06,0.04,0.02])

        # define perturbation coefficient for calculating Jacobian matrix
        # the direction is determine by random number generator either 1 or -1
        # this coefficient will be evaluated based on the cost function chi2
        # if chi2 reduces, this coefficient should be reduced too
        if i == 0:
            p_tau = c_ptau[0]*[-1,1][random.randrange(2)]*xi_arr[i,0]   
            p_reff = c_ptau[0]*[-1,1][random.randrange(2)]*xi_arr[i,1]
        
        # print statement
        print("Info     | Perturbation coefficient p_tau = %.2f and p_reff = %.2f" %(p_tau,p_reff))
        
        # define cloud filename
        cloud_file1 = dummy_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        cloud_file2 = dummy_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        cloud_file3 = dummy_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)
            
        # define lwc :: in the end it would be modified as we set tau550
        lwc = 0.2
        
        # define cloud profiles
        cloud_profile(cloud_file1,float(ztop),float(zbase),lwc,xi_arr[i,1])
        cloud_profile(cloud_file2,float(ztop),float(zbase),lwc,xi_arr[i,1])
        cloud_profile(cloud_file3,float(ztop),float(zbase),lwc,xi_arr[i,1]+p_reff)
        
        # +++++++++++++
        # forward model
        # +++++++++++++
        # input-output filenames
        inp_file1 = dummy_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        inp_file2 = dummy_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        inp_file3 = dummy_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)

        out_file1 = dummy_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        out_file2 = dummy_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        out_file3 = dummy_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)
        
        # generate input
        uvspec_input(lib_data_dir=lib_data_dir,cloud_file=cloud_file1,input_file=inp_file1,\
                     sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,albedo_file=albedo_file,lambda0=lambda0_fwd,lambda1=lambda1_fwd,\
                     tau550=str(xi_arr[i,0]),aerosol_season=aerosol_season,aerosol_haze=aerosol_haze)

        uvspec_input(lib_data_dir=lib_data_dir,cloud_file=cloud_file2,input_file=inp_file2,\
                     sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,albedo_file=albedo_file,lambda0=lambda0_fwd,lambda1=lambda1_fwd,\
                     tau550=str(xi_arr[i,0]+p_tau),aerosol_season=aerosol_season,aerosol_haze=aerosol_haze)

        uvspec_input(lib_data_dir=lib_data_dir,cloud_file=cloud_file3,input_file=inp_file3,\
                     sza=sza,phi0=phi0,phi=phi,zout=zout,doy=doy,albedo_file=albedo_file,lambda0=lambda0_fwd,lambda1=lambda1_fwd,\
                     tau550=str(xi_arr[i,0]),aerosol_season=aerosol_season,aerosol_haze=aerosol_haze)

        # define run script
        run_script1 = "uvspec <" + inp_file1 + "> " + out_file1
        run_script2 = "uvspec <" + inp_file2 + "> " + out_file2
        run_script3 = "uvspec <" + inp_file3 + "> " + out_file3
        
        # run libradtran
        os.system(run_script1)
        os.system(run_script2)
        os.system(run_script3)

        # +++++++++++++++++++++
        # reading uvspec output
        # +++++++++++++++++++++
        # define lbl outfile
        out_lbl_file1 = dummy_dir + "lbl_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        out_lbl_file2 = dummy_dir + "lbl_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])         
        out_lbl_file3 = dummy_dir + "lbl_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)         
        
        # read output 
        read_uvspec_output(input_file=out_file1,output_file=out_lbl_file1)
        read_uvspec_output(input_file=out_file2,output_file=out_lbl_file2)
        read_uvspec_output(input_file=out_file3,output_file=out_lbl_file3)
        
        # ++++++++++++++++++++
        # spectral convolution
        # ++++++++++++++++++++
        # define output filename
        out_conv_file1 = dummy_dir + "spectra_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        out_conv_file2 = dummy_dir + "spectra_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])        
        out_conv_file3 = dummy_dir + "spectra_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)        

        # convolution        
        l1b_convolution(input_file=out_lbl_file1,spectral_grid=lambda_meas,output_file=out_conv_file1)
        l1b_convolution(input_file=out_lbl_file2,spectral_grid=lambda_meas,output_file=out_conv_file2)
        l1b_convolution(input_file=out_lbl_file3,spectral_grid=lambda_meas,output_file=out_conv_file3)

        # ++++++++++++++++++
        # calculate jacobian
        # ++++++++++++++++++
        # define outfile
        out_kmat_file = out_dir + "kmat_" + str(pix_id).zfill(5) + "_" + str(i+1).zfill(2) + ".dat"
        
        # call kmat module
        kmat_arr = kmat_module(input_file1=out_conv_file1,input_file2=out_conv_file2,input_file3=out_conv_file3,\
                               tau1=xi_arr[i,0],reff1=xi_arr[i,1],tau2=xi_arr[i,0]+p_tau,reff2=xi_arr[i,1]+p_reff,\
                               lambda_meas=lambda_meas,rad_meas=rad_meas,sy=np.diag(sy),outfile=out_kmat_file)

        # define modeled spectra
        rad_fwd = kmat_arr[:,3]

        # define kmat array
        kmat = kmat_arr[:,4:]
        
        # +++++++++++++++++++++++++++++++++++++++++++++++
        # optimal estimation module :: see Rodgers (2000)
        # +++++++++++++++++++++++++++++++++++++++++++++++
        # for i == 0 :: starting point
        if i == 0:
            # regularization parameter :: initial value
            gamma = np.mean(np.diag(np.matmul(kmat.T,np.matmul(inv(sy),kmat))))
            
            # define cost function :: initial value
            chi2_0 = np.matmul((rad_meas-rad_fwd).T,np.matmul(inv(sy),rad_meas-rad_fwd)) + np.matmul((x0-xa).T,np.matmul(inv(sa),x0-xa))

            # assign chi2 to array
            chi2_arr[i] = chi2_0                       

            # retrieval module
            (xi,sx,G,A,chi2) = l2m_retrieval(rad_meas=rad_meas,rad_fwd=rad_fwd,sy=sy,sa=sa,kmat=kmat,xii=xi_arr[i,:],xa=xa,gamma=gamma)
    
            # assign chi2
            chi2_arr[i+1] = chi2
        
        # for general condition
        elif i > 0:
            # retrieval module
            (xi,sx,G,A,chi2) = l2m_retrieval(rad_meas=rad_meas,rad_fwd=rad_fwd,sy=sy,sa=sa,kmat=kmat,xii=xi_arr[i,:],xa=xa,gamma=gamma)
    
            # assign chi2
            chi2_arr[i+1] = chi2

        # ++++++++++++++++++++++++++
        # evaluation of current step
        # ++++++++++++++++++++++++++
        # coefficient to update gamma :: e.g., Rodgers (2000)
        f = 10          
        
        # initial condition
        if i == 0:
            # update state vector
            xi_arr[i+1,:] = xi

        # general condition            
        elif i > 0 and chi2_arr[i+1] > chi2_arr[i]:
            # define new gamma
            gamma = gamma*f     # gradient descent
            
            # define perturbation coefficient
            p_tau = c_ptau[0]*[-1,1][random.randrange(2)]*xi_arr[i,0]   
            p_reff = c_ptau[0]*[-1,1][random.randrange(2)]*xi_arr[i,1]
            
            # update state vector using Gain matrix :: linear approximation :: shortcut
            dx = np.matmul(G,(rad_meas-rad_fwd))

            # updata state vector            
            xi_arr[i+1,:] = xi_arr[i,:] + dx

        # semi trusted region            
        elif i > 0 and chi2_arr[i+1] > nmeas and chi2_arr[i+1] < chi2_arr[i]:
            # define new gamma
            gamma = gamma/f     # gradient descent

            # define delta state for the reference when calculating Jacobian matrix
            p_tau = c_ptau[1]*[-1,1][random.randrange(2)]*xi_arr[i,0]   
            p_reff = c_ptau[1]*[-1,1][random.randrange(2)]*xi_arr[i,1]
            
            # update state vector using Gain matrix :: linear approximation :: shortcut
            dx = np.matmul(G,(rad_meas-rad_fwd))
            
            # update state vector
            xi_arr[i+1,:] = xi + dx

        # trusted region        
        elif i > 0 and chi2_arr[i+1] <= nmeas:
            # define new gamma
            gamma = 0           # pure Gauss-Newton approach

            # define delta state for the reference when calculating Jacobian matrix
            p_tau = c_ptau[1]*[-1,1][random.randrange(2)]*xi_arr[i,0]   
            p_reff = c_ptau[1]*[-1,1][random.randrange(2)]*xi_arr[i,1]
            
            # update state vector using Gain matrix :: linear approximation
            dx = np.matmul(G,(rad_meas-rad_fwd))
            
            # update state vector
            xi_arr[i+1,:] = xi + dx            
        
        # print statement
        print("Info     | tau = %.3f and reff = %.3f" %(xi_arr[i+1,0],xi_arr[i+1,1]))
        print("Info     | cost function chi2 = %.3f" %(chi2_arr[i+1])) 
    
        # ++++++++++++++++++++
        # retrieval diagnostic
        # ++++++++++++++++++++
        # define messages:        
        # 1     = converged
        # 2     = does not converged :: exceed maximum number of iteration
        # 3     = does not converged :: boundary hit
        
        # define ratio chi2
        r_chi2 = np.abs(1 - (chi2_arr[i+1]/chi2_arr[i]))
        
        # define threshold (ratio in percent) between chi2_i+1 and chi2_i+1
        # a retrieval is converged when the change on chi2 is neggligible
        # don't use absolute threshold as it changes with n_measurements
        t_dchi2 = 0.015  
        
        # 1 = converged :: should be more than one iteration
        if i > 0 and r_chi2 < t_dchi2:     
            # print statement
            print("Info     | Convergence ID = 1 :: Retrieval succeed!")
            
            # define filenames
            stat_file = out_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = out_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=1,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)
        
            # stop iteration
            break
            
        # 2 = does not converged :: exceed maximum number of iteration        
        elif i == int(iter_max)-1:      
            # print statement
            print("Info     | Convergence ID = 2 :: Exceeding number maximum iteration!")

            # define filenames
            stat_file = out_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = out_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=2,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)
        
            # stop iteration
            break
            
        # 3 = does not converged :: boundary hit
        elif xi_arr[i+1,0] <= tau_min or xi_arr[i+1,0] >= tau_max or \
        xi_arr[i+1,1] <= reff_min or xi_arr[i+1,1] >= reff_max or \
        xi_arr[i,0]+p_tau <= tau_min or xi_arr[i,0]+p_tau >= tau_max or \
        xi_arr[i,1]+p_reff <= reff_min or xi_arr[i,1]+p_reff >= reff_max or \
        xi_arr[i,0] <= tau_min or xi_arr[i,0] >= tau_max or \
        xi_arr[i,1] <= reff_min or xi_arr[i,1] >= reff_max:
            # print statement
            print("Info     | Convergence ID = 3 :: Lower or upper boundary hit!")

            # define filenames
            stat_file = out_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = out_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=3,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)

            # stop iteration
            break
            
        # update iteration
        i += 1

