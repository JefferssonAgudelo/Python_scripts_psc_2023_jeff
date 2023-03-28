#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 10:44:43 2023

@author: jeffersson_agudelo
"""


import glob, os
import numpy as np
import matplotlib.pyplot as plt
import h5py
#import vaex
import gc
gc.enable()


# load the file
############################################################################################################################
# To work with fields outputs 
############################################################################################################################
# set the path to where the data are
path = '/Volumes/PSC_CORI_DATA/TESTS_2022/New_V_tests/runs_with_L_05/testing_dB/psc_harris_dummy_2'
os.chdir(path)

# Fields file
file_f=sorted(glob.glob(os.path.join(path, 'pfd.005650_p000000.h5'))) # This loads all the files in the directory with the .h5 extension 
# Moments file
file_m=sorted(glob.glob(os.path.join(path, 'pfd_moments.005650_p000000.h5'))) # This loads all the files in the directory with the .h5 extension 

a=0; 

# Load the fields variables e,h,j
###############################################################################
#------------------------------------------------------------------------------
filename1_f= file_f[a]
print(filename1_f)
h5_1st_f = h5py.File(filename1_f, 'r') # Load the file but not the data yet.
list(h5_1st_f.keys()) # This helps you to know what index you want [15]

#------------------------------------------------------------------------------
jeh_fields = h5_1st_f[list(h5_1st_f.keys())[15]] ; #Like this the outcome is a group
#------------------------------------------------------------------------------

# Electric field
#------------------------------------------------------------------------------ 
ex_ec = jeh_fields['ex_ec']; ex_ec=ex_ec['p0']; ex_ec=ex_ec['3d']; ex_ec=ex_ec[:,:,:]; ex_ec=np.transpose(ex_ec, (2, 1, 0)) # This (z,y,x) -> (x,y,z)
ey_ec = jeh_fields['ey_ec']; ey_ec=ey_ec['p0']; ey_ec=ey_ec['3d']; ey_ec=ey_ec[:,:,:]; ey_ec=np.transpose(ey_ec, (2, 1, 0))
ez_ec = jeh_fields['ez_ec']; ez_ec=ez_ec['p0']; ez_ec=ez_ec['3d']; ez_ec=ez_ec[:,:,:]; ey_ec=np.transpose(ey_ec, (2, 1, 0))
#------------------------------------------------------------------------------

# Magnetic field
#------------------------------------------------------------------------------ 
hx_fc = jeh_fields['hx_fc']; hx_fc=hx_fc['p0']; hx_fc=hx_fc['3d']; hx_fc=hx_fc[:,:,:]; hx_fc=np.transpose(hx_fc, (2, 1, 0)) # This (z,y,x) -> (x,y,z)
hy_fc = jeh_fields['hy_fc']; hy_fc=hy_fc['p0']; hy_fc=hy_fc['3d']; hy_fc=hy_fc[:,:,:]; hy_fc=np.transpose(hy_fc, (2, 1, 0))
hz_fc = jeh_fields['hz_fc']; hz_fc=hz_fc['p0']; hz_fc=hz_fc['3d']; hz_fc=hz_fc[:,:,:]; hy_fc=np.transpose(hy_fc, (2, 1, 0))
#------------------------------------------------------------------------------

# Electric current field
#------------------------------------------------------------------------------
jx_ec = jeh_fields['jx_ec']; jx_ec=jx_ec['p0']; jx_ec=jx_ec['3d']; jx_ec=jx_ec[:,:,:]; jx_ec=np.transpose(jx_ec, (2, 1, 0)) # This (z,y,x) -> (x,y,z)
jy_ec = jeh_fields['jy_ec']; jy_ec=jy_ec['p0']; jy_ec=jy_ec['3d']; jy_ec=jy_ec[:,:,:]; jy_ec=np.transpose(jy_ec, (2, 1, 0))
jz_ec = jeh_fields['jz_ec']; jz_ec=jz_ec['p0']; jz_ec=jz_ec['3d']; jz_ec=jz_ec[:,:,:]; jy_ec=np.transpose(jy_ec, (2, 1, 0))
#------------------------------------------------------------------------------
###############################################################################


# Load the moment variables 
###############################################################################
#------------------------------------------------------------------------------
filename1_m= file_m[a]
print(filename1_m)
h5_1st_m = h5py.File(filename1_m, 'r') # Load the file but not the data yet.
list(h5_1st_m.keys()) # This helps you to know what index you want 


#------------------------------------------------------------------------------
all_1st = h5_1st_m[list(h5_1st_m.keys())[0]] ; #Like this the outcome is a group
#------------------------------------------------------------------------------

# Run this to see the name of the variables that you want to load
#------------------------------------------------------------------------------
list(all_1st)
#------------------------------------------------------------------------------
# In my case my especies are i_UP, i_BA, e_UP, e_BA. Change yours accordingly 

#Momentum (p = gamma * m_s * v) components ions UP# 
#------------------------------------------------------------------------------
px_i_UP=all_1st['px_i_UP'] ; px_i_UP=px_i_UP['p0'] ; px_i_UP=px_i_UP['3d'] ; px_i_UP=px_i_UP[:,:,:] ; px_i_UP=np.transpose(px_i_UP, (2, 1, 0))
py_i_UP=all_1st['py_i_UP'] ; py_i_UP=py_i_UP['p0'] ; py_i_UP=py_i_UP['3d'] ; py_i_UP=py_i_UP[:,:,:] ; py_i_UP=np.transpose(py_i_UP, (2, 1, 0))
pz_i_UP=all_1st['pz_i_UP'] ; pz_i_UP=pz_i_UP['p0'] ; pz_i_UP=pz_i_UP['3d'] ; pz_i_UP=pz_i_UP[:,:,:] ; pz_i_UP=np.transpose(pz_i_UP, (2, 1, 0))
#------------------------------------------------------------------------------

#Momentum (p = gamma * m_s * v) components electrons UP# 
#------------------------------------------------------------------------------
px_e_UP=all_1st['px_e_UP'] ; px_e_UP=px_e_UP['p0'] ; px_e_UP=px_e_UP['3d'] ; px_e_UP=px_e_UP[:,:,:] ; px_e_UP=np.transpose(px_e_UP, (2, 1, 0))
py_e_UP=all_1st['py_e_UP'] ; py_e_UP=py_e_UP['p0'] ; py_e_UP=py_e_UP['3d'] ; py_e_UP=py_e_UP[:,:,:] ; py_e_UP=np.transpose(py_e_UP, (2, 1, 0))
pz_e_UP=all_1st['pz_e_UP'] ; pz_e_UP=pz_e_UP['p0'] ; pz_e_UP=pz_e_UP['3d'] ; pz_e_UP=pz_e_UP[:,:,:] ; pz_e_UP=np.transpose(pz_e_UP, (2, 1, 0))
#------------------------------------------------------------------------------

#Charge density (q_{s}*n_{s}) ions and electrons UP
#------------------------------------------------------------------------------
n_i_UP=all_1st['rho_i_UP'] ; n_i_UP=n_i_UP['p0'] ; n_i_UP=n_i_UP['3d'] ; n_i_UP=n_i_UP[:,:,:] ; n_i_UP=np.transpose(n_i_UP, (2, 1, 0))
n_e_UP=all_1st['rho_e_UP'] ; n_e_UP=n_e_UP['p0'] ; n_e_UP=n_e_UP['3d'] ; n_e_UP=n_e_UP[:,:,:] ; n_e_UP=np.transpose(n_e_UP, (2, 1, 0))
#------------------------------------------------------------------------------
###############################################################################



#Make the operations that you need
#------------------------------------------------------------------------------

# Electromagnetic energy
#------------------------------------------------------------------------------
U_em = 0.5 * (pow(ex_ec,2) + pow(ey_ec,2) + pow(ez_ec,2)) + \
    0.5 * (pow(hx_fc,2) + pow(hy_fc,2) + pow(hz_fc,2)) 
#------------------------------------------------------------------------------


# take a 2D slice
#------------------------------------------------------------------------------
U_em_yz = U_em[0,:,:]
#------------------------------------------------------------------------------


#Make plots
###############################################################################
fig, ax0 = plt.subplots()
plt.pcolor(U_em_yz)
ax0.set_title('U em')
#Make plots
###############################################################################

