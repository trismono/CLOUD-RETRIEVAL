#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program for

Created on %(date)

@author : trismonock
@mail   : T.C.Krisna@sron.nl
"""

# import standard modules
import os, numpy as np, random, time
from numpy.linalg import inv

# import local routines
from read_config import read_ret_config, read_geometry_config
from uvspec_input import uvspec_input
from cloud_profile import cloud_profile
from convolution import convolution
from read_uvspec_output import read_uvspec_output
from write_spectra_to_ascii import write_spectra_to_ascii
from kmat_module import kmat_module
from l2m_retrieval import l2m_retrieval
from diagnostic_module import diagnostic

def l2m_wrapper(base_dir, libradtran_dir, pix_id):     
    # ++++++++++++
    # define paths
    # ++++++++++++
    input_dir = os.path.realpath(base_dir + "../INPUT") + "/" 
    output_dir = os.path.realpath(base_dir + "../OUTPUT") + "/"
    case_dir = os.path.realpath(base_dir + "../../01_CREATE_CASE/OUTPUT") + "/"
    l1b_dir = os.path.realpath(base_dir + "../../03_L1B/OUTPUT") + "/"
    lib_data_dir = os.path.realpath(libradtran_dir + "data")               # libtradtran internal database
    
    # +++++++++++++++++++++
    # read retrieval config
    # +++++++++++++++++++++
    # define config file
    ret_config_file = input_dir + "l2m_config.dat"
    
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
    
    # read radiance
    if flag_noise == "0":       # without noise
        rad_meas_vnir = meas_vnir_data[:,2]
        rad_meas_swir = meas_swir_data[:,2] 
    elif flag_noise == "1":     # with noise
        rad_meas_vnir = meas_vnir_data[:,1]
        rad_meas_swir = meas_swir_data[:,1]
 

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
    nstate = 2       # number of state to be retrieved: tau and reff
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
    xi_arr = np.zeros(shape=(int(iter_max)+1,nstate))
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

    # define time start
    time0 = time.perf_counter()
    
    # start iteration : the number of iteration is set in the retrieval setting
    i = 0
    while i < int(iter_max):
        # print statement
        print("Info         | Pixel ID %i :: Iteration %i " %(pix_id,i+1))

        # +++++++++++++++++++++        
        # running forward model
        # +++++++++++++++++++++
        # define perturbation coefficient in percent
        c_ptau = np.array([0.04,0.02,0.01])
        c_preff = np.array([0.04,0.02,0.01])

        # define perturbation coefficient for calculating Jacobian matrix
        # the direction is determine by random number generator either 1 or -1
        # this coefficient will be evaluated based on the cost function chi2
        # if chi2 reduces, this coefficient should be reduced too
        if i == 0:
            p_tau = c_ptau[0]*[-1,1][random.randrange(2)]*x0[0]   
            p_reff = c_preff[0]*[-1,1][random.randrange(2)]*x0[1]
        
        # print statement
        print("Info         | Pixel ID %i :: Perturbation coefficient p_tau = %.2f and p_reff = %.2f micron" %(pix_id,p_tau,p_reff))
        
        # define cloud filename
        cloud_file1 = input_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        cloud_file2 = input_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        cloud_file3 = input_dir + "cloud_%s_%.3f_%.3f.dat" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)
            
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
        inp_file1 = input_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        inp_file2 = input_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        inp_file3 = input_dir + "cloud_%s_%.3f_%.3f.inp" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)

        out_file1 = input_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1])
        out_file2 = input_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0]+p_tau,xi_arr[i,1])
        out_file3 = input_dir + "cloud_%s_%.3f_%.3f.out" %(str(pix_id).zfill(5),xi_arr[i,0],xi_arr[i,1]+p_reff)
        
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
        # read output 
        data_lbl1 = read_uvspec_output(input_file=out_file1)
        data_lbl2 = read_uvspec_output(input_file=out_file2)
        data_lbl3 = read_uvspec_output(input_file=out_file3)

        # define lbl spectra outfile
        out_lbl_file1 = output_dir + "lbl_spectra_%s_%s.dat" %(str(pix_id).zfill(5),str(i+1).zfill(2))
                
        # write lbl spectra
        write_spectra_to_ascii(filename=out_lbl_file1, data=data_lbl1)
        
        # ++++++++++++++++++++
        # spectral convolution
        # ++++++++++++++++++++
        # convolution        
        data_conv1 = convolution(data_fwd=data_lbl1, spectral_grid=lambda_meas)
        data_conv2 = convolution(data_fwd=data_lbl2, spectral_grid=lambda_meas)
        data_conv3 = convolution(data_fwd=data_lbl3, spectral_grid=lambda_meas)

        # define output filename
        out_conv_file1 = output_dir + "spectra_%s_%s.dat" %(str(pix_id).zfill(5),str(i+1).zfill(2))

        # write convolved spectra
        write_spectra_to_ascii(filename=out_conv_file1, data=data_conv1)
        
        # ++++++++++++++++++
        # calculate jacobian
        # ++++++++++++++++++
        # define outfile
        out_kmat_file = output_dir + "kmat_" + str(pix_id).zfill(5) + "_" + str(i+1).zfill(2) + ".dat"
        
        # call kmat module
        kmat_arr = kmat_module(data_ref=data_conv1,data_ptau=data_conv2,data_preff=data_conv3,\
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
            # regularization parameter :: initial value - will be evaluated in each iteration
            gamma = 1000
            
            # define cost function :: initial value
            chi2_0 = np.matmul((rad_meas-rad_fwd).T,np.matmul(inv(sy),rad_meas-rad_fwd)) + np.matmul((x0-xa).T,np.matmul(inv(sa),x0-xa))

            # assign chi2 to array
            chi2_arr[i] = chi2_0                       

            # retrieval module
            (xi,sx,G,A,chi2) = l2m_retrieval(rad_meas=rad_meas,rad_fwd=rad_fwd,sy=sy,sa=sa,kmat=kmat,xi=xi_arr[i,:],xa=xa,gamma=gamma)
    
            # assign chi2
            chi2_arr[i+1] = chi2
        
        # for general condition
        elif i > 0:
            # retrieval module
            (xi,sx,G,A,chi2) = l2m_retrieval(rad_meas=rad_meas,rad_fwd=rad_fwd,sy=sy,sa=sa,kmat=kmat,xi=xi_arr[i,:],xa=xa,gamma=gamma)
    
            # assign chi2
            chi2_arr[i+1] = chi2

        # ++++++++++++++++++++++++++
        # evaluation of current step
        # ++++++++++++++++++++++++++
        # coefficient to update gamma 
        f = 3          
        
        # initial condition
        if i == 0:
            # define perturbation coefficient
            p_tau = c_ptau[0]*[-1,1][random.randrange(2)]*xi[0]   
            p_reff = c_preff[0]*[-1,1][random.randrange(2)]*xi[1]
            
            # update state vector
            xi_arr[i+1,:] = xi

        # general condition            
        elif i > 0 and chi2_arr[i+1] > chi2_arr[i]:
            # define new gamma
            gamma = gamma*f     # gradient descent
            
            # define perturbation coefficient
            p_tau = c_ptau[0]*[-1,1][random.randrange(2)]*xi[0]   
            p_reff = c_preff[0]*[-1,1][random.randrange(2)]*xi[1]

            # updata state vector            
            xi_arr[i+1,:] = xi_arr[i,:]

        # semi trusted region            
        elif i > 0 and chi2_arr[i+1] > nmeas and chi2_arr[i+1] < chi2_arr[i]:
            # define new gamma
            gamma = gamma/f     # gradient descent

            # define delta state for the reference when calculating Jacobian matrix
            p_tau = c_ptau[1]*[-1,1][random.randrange(2)]*xi[0]  
            p_reff = c_preff[1]*[-1,1][random.randrange(2)]*xi[1]
            
            # update state vector
            xi_arr[i+1,:] = xi

        # trusted region        
        elif i > 0 and chi2_arr[i+1] <= nmeas:
            # define new gamma
            gamma = 1           # pure Gauss-Newton approach

            # define delta state for the reference when calculating Jacobian matrix
            p_tau = c_ptau[2]*[-1,1][random.randrange(2)]*xi[0]   
            p_reff = c_preff[2]*[-1,1][random.randrange(2)]*xi[1]
            
            # update state vector
            xi_arr[i+1,:] = xi            
        
        # print statement
        print("Info         | Pixel ID %i :: State tau = %.3f and reff = %.3f micron" %(pix_id,xi_arr[i+1,0],xi_arr[i+1,1]))
        print("Info         | Pixel ID %i :: Cost function chi2 = %.3f" %(pix_id,chi2_arr[i+1])) 
    
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
        t_dchi2 = 0.05
        
        # 1 = converged :: should be more than one iteration
        if i > 0 and r_chi2 < t_dchi2 and chi2_arr[i+1] <= nmeas:     
            # print statement
            print("Info         | Pixel ID %i :: Convergence ID = 1 Retrieval succeed!" %(pix_id))
            
            # define filenames
            stat_file = output_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = output_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=1,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)
        
            # stop iteration
            break
            
        # 2 = does not converged :: exceed maximum number of iteration        
        elif i == int(iter_max)-1:      
            # print statement
            print("Info         | Pixel ID %i :: Convergence ID = 2 Exceeding number maximum iteration!" %(pix_id))

            # define filenames
            stat_file = output_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = output_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=2,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)
        
            # stop iteration
            break
            
        # 3 = does not converged :: boundary hit
        elif xi[0] <= tau_min or xi[0] >= tau_max or \
        xi[1] <= reff_min or xi[1] >= reff_max or \
        xi[0]+p_tau <= tau_min or xi[0]+p_tau >= tau_max or \
        xi[1]+p_reff <= reff_min or xi[1]+p_reff >= reff_max:
            # print statement
            print("Info         | Pixel ID %i :: Convergence ID = 3 Lower or upper boundary hit!" %(pix_id))

            # define filenames
            stat_file = output_dir + "statistics_" + str(pix_id).zfill(5) + ".dat"
            gain_file = output_dir + "gain_" + str(pix_id).zfill(5) + ".dat"
            
            # writing statistics file
            diagnostic(stat_file=stat_file,pix_id=pix_id,message=3,niter=i+1,chi2=chi2_arr[i+1],x=xi_arr[i+1,:],gain_file=gain_file,lambda_meas=lambda_meas,G=G)

            # stop iteration
            break
            
        # update iteration
        i += 1

    # define end time
    time1 = time.perf_counter()
    
    # calculate execution time
    time_total = time1-time0
    time_mean = time_total/(i+1)
    
    # print statement
    print("Info         | Pixel ID %i :: Retrieval done in %.2f sec" %(pix_id,time_total))
    print("Info         | Pixel ID %i :: Elapsed time per iteration = %.2f sec" %(pix_id,time_mean))
