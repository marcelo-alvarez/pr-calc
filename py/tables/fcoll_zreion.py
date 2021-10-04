#!/usr/bin/env python
# coding: utf-8

# Imports

import sys
#import Pk
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import peaks
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.integrate import quad
from colossus.lss import mass_function
import matplotlib.pyplot as plt
import pickle

with open('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/parameterfiles/param.col', 'rb') as handle:
    b = pickle.load(handle)

print(b)
print(b["flat"])
# Variables + cosmology

#we set our own cosmology
params         = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.96}
params         = {'flat': b["flat"], 'H0': float(b["H0"]), 'Om0': float(b["Om0"]), 'Ob0': float(b["Ob0"]), 'sigma8': float(b["sigma8"]), 'ns': 0.96}
cosmo          = cosmology.setCosmology('myCosmo', params)
h              = cosmo.H0 / 100.

# Mmin and zeta

Mmin_Msun      = float(sys.argv[1]) #Mmin in Msun
print(Mmin_Msun, 'Mmin value as used in the python code')
Mmax           = 1E15     
zeta           = float(sys.argv[2])

# create grid of delta_R and R and the zreion table
# can mess with the grid later but an initial guess follows
nR             = 100
ndelta_R       = 101 # odd ndelta_R for a delta_R=0 entry with symmetric range
Rvalsmin       = 1.0
Rvalsmax       = 1e3
delta_Rmin     = -1
delta_Rmax     =  1
#R should be in Msun/h. I'm not sure we have done that
Rvals          = np.logspace(np.log10(Rvalsmin), np.log10(Rvalsmax), nR)
delta_Rvals    = np.linspace(delta_Rmin, delta_Rmax, ndelta_R)
zreion_table   = np.zeros((nR, ndelta_R))

# create grid of z and the fcoll table
dz             = 0.5   
zmin           = 0.0 # just for completeness
zmax           = 20  # most reionization models we don't care about z>20 
nz             = int((zmax - zmin) / dz)
fcoll_table    = np.zeros(nz)

rho_Mpc        = 2.775E11 * cosmo.Om0 * h**2
 
zvals          = np.linspace(zmin, zmax, nz)  # zvals will be approximately spaced at dz
fcoll_table    = np.zeros((nz, nR, ndelta_R)) # holds all the collapsed fractions

uncond_fcoll   = np.zeros(len(zvals))
uncond_fcollPS = np.zeros(len(zvals))

# meshgrids for independent variables
Rvals2d, delta_Rvals2d          = np.meshgrid(Rvals,   delta_Rvals, indexing='ij')
zvals3d, Rvals3d, delta_Rvals3d = np.meshgrid(zvals,   Rvals,   delta_Rvals, indexing='ij')

# convert Mmin to Rmin
Rmin_Mpc       = (3*Mmin_Msun/4/np.pi/rho_Mpc)**(1./3.)   #R=(3M/4 pi)**1/3
Rmin           = Rmin_Mpc*h # Rmin is used by colussus and is in Msun/h

dRzero_index   = int((ndelta_R+1)/2)-1

# hmf and halo mass definition
hmf_mdef       = 'fof'
hmf_model      = 'press74'
#hmf_mdef      = 'fof'
#hmf_model     = 'sheth99'
#hmf_mdef      = '200m'
#hmf_model     = 'tinker08'

# plotting 
plt.rcParams['figure.figsize'] = [12, 8]


# Function to compute PS conditional fcoll

# function for fcoll = PS conditional collapsed fraction (ccfPS)
def ccfPS(z, R, delta_R):
    delta_c   = peaks.collapseOverdensity(corrections = True,z = z)
    sigma_R   = cosmo.sigma(R, z)  #R in Mpc/h. MI I'm not sure R is in Msun/h yet
    sigma_min = cosmo.sigma(Rmin, z) #Rmin in Mpc/h
    arg       = ((delta_c - delta_R) / (np.sqrt(2.*sigma_min**2 - sigma_R**2)))
    f         = erfc(arg)
    return f
ccfPS = np.vectorize(ccfPS)


# Function to compute unconditional fcoll from colossus for a given hmf model and mass definition

def integrand(lnM, rho_Mpc, redshift, model, mdef):
    '''
    This function computes the collapse fraction as an integral of the mass function indicated by the user
    '''
    h  = cosmo.H0 / 100. 
    M  = np.exp(lnM)     # input mass assumed to be in Msun
    M /= h               # convert from Msun to Msun/h for Colussus
    
    dndlnM = mass_function.massFunction(M, redshift, mdef = mdef, model = model, q_out = 'dndlnM')

    dndlnM *= h**3 # convert mass function from 1/(Mpc/h)^3 to 1/Msun/Mpc^3
    M      *= h    # convert mass back to Msun from Msun/h 
    
    dfcolldlnM = 1/rho_Mpc * M * dndlnM   # convert dfcoll/dlnM = 1/rho*M*dn/dlnM

    return dfcolldlnM

def ucf(z,hmf_model,hmf_mdef):
    # integrate hmf to obtain collapsed fraction 
    from scipy import integrate
    uncond_fcoll,   err = quad(integrand, np.log(Mmin_Msun), np.log(Mmax), 
                                args=(rho_Mpc, z, hmf_model, hmf_mdef))
    return uncond_fcoll
ucf = np.vectorize(ucf)


# Function to compute conditional fcoll using equation 3 from overleaf

def ccf(z, R, delta_R, ucf_func, ucfPS_func):
    # eqn 3 from overleaf 
    # cond_fcoll(z,R,delta_R) = uncond_fcoll(z) / uncond_fcollPS * cond_fcollPS
    return ucf_func(z) / ucfPS_func(z) * ccfPS(z, R, delta_R) 


# Function to compute reionization redshift at a given radius and density contrast

def zreion(R, delta_R, ucf_func, ucfPS_func, zeta):
    # note zvals is a global variable!!!
    fcoll_zvals = ccf(zvals,R,delta_R,ucf_func, ucfPS_func)
    zeta_times_fcoll_zvals = zeta * fcoll_zvals
    z_of_zeta_times_fcoll = interp1d(zeta_times_fcoll_zvals, zvals, fill_value="extrapolate")
    zreion = z_of_zeta_times_fcoll(1.0)
    return zreion
zreion = np.vectorize(zreion)


# Populate collapsed fraction tables and set up collapsed fraction interpolation tables
# tables for fcoll(z)
#   these are the unconditional collapsed fractions for PS 
#   as well as that for the specified hmf model
uncond_fcoll   = ucf(zvals,hmf_model,hmf_mdef).astype(np.float32)
uncond_fcollPS = ucf(zvals,'press74','fof').astype(np.float32)
ucf_interp     = interp1d(zvals, uncond_fcoll)
ucfPS_interp   = interp1d(zvals, uncond_fcollPS)


# Populate zreion table

# table for zreion(R,delta_R)
#   this is the reionization redshift at a given overdensity and scale 
#   from the condition
#   zeta * fcoll(R,delta_R,zreion) = 1
zreion_table = zreion(Rvals2d, delta_Rvals2d, ucf_interp, ucfPS_interp, zeta).astype(np.float32)

# Write tables in binaries with minimal headers easily readable in C

# fcoll(z) (unconditional mass function)  ......why writing the unconditional to file?
fcoll_tablefile='fcoll_table.tab'
f = open(fcoll_tablefile,'wb')
zrange       = np.asarray([zvals[0],zvals[-1]]).astype(np.float32)
dim          = np.asarray([len(zvals)]).astype(np.int32)
fcoll_table  = uncond_fcoll.astype(np.float32)
zrange.tofile(f)
dim.tofile(f)
fcoll_table.tofile(f)
f.close()

# zreion(R,delta_R) (reionization redshifts) 
zreion_tablefile='zreion_table.tab'
f = open(zreion_tablefile,'wb')
ranges = [Rvalsmin,Rvalsmax,delta_Rmin,delta_Rmax]
ranges = np.asarray(ranges).astype(np.float32)
dims   = [nR, ndelta_R]
dims   = np.asarray(dims).astype(np.int32)
ranges.tofile(f)
dims.tofile(f)
zreion_table.tofile(f)
f.close()


# Done with table generation. Functions below are for reading tables and testing results.

# function to read zreion table
def read_zreiontable(fname):
    f = open(fname,'rb')  #fname has to be defined somehow so it knows what file to open, right?
    Rrange       = np.fromfile(f,dtype=np.float32, count=2)
    deltaRrange  = np.fromfile(f,dtype=np.float32, count=2)
    dims         = np.fromfile(f,dtype=np.int32,   count=2)
    zreion_table = np.fromfile(f,dtype=np.float32, count=dims[0]*dims[1])
    zreion_table = np.reshape(zreion_table,(dims[0],dims[1]))
    print('read table of dimensions',dims)
    return Rrange, deltaRrange, dims, zreion_table

# check table read and zreion
def plot_zreion_from_table(fname):
    #plot z_reion 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
    Rrange, deltaRrange, dims, zreion_table = read_zreiontable(fname)
    Rvals      = np.logspace(np.log10(Rrange[0]),     np.log10(Rrange[1]),     dims[0])
    deltaRvals = np.logspace(         deltaRrange[0],          deltaRrange[1], dims[1])
    plt.plot(Rvals, zreion_table[:,0]           ,'r:',label='$\delta_R=$'+str(delta_Rvals[0]))
    plt.plot(Rvals, zreion_table[:,dRzero_index],'k-',label='$\delta_R=$'+str(delta_Rvals[dRzero_index]))
    plt.plot(Rvals, zreion_table[:,-1]          ,'b:',label='$\delta_R=$'+str(delta_Rvals[-1]))

    plt.legend(fontsize=20)
    plt.tick_params(labelsize=20)
    plt.xscale('log')
    plt.xlabel('R [Mpc]',fontsize=20)
    plt.ylabel(r'$z_{\rm reion}$',fontsize=20)
    plt.savefig('zreion.png',bbox_inches='tight')
    #plt.show()
    print('zreion should increase towards higher R and delta_R')
    print('this is not an exhaustive test, just a single quick check')
    
def show_zreion_from_table(fname):
    #plot z_reion 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
    Rrange, deltaRrange, dims, zreion_table = read_zreiontable(fname)
    Rvals      = np.logspace(np.log10(Rrange[0]),     np.log10(Rrange[1]),     dims[0])
    deltaRvals = np.logspace(         deltaRrange[0],          deltaRrange[1], dims[1])
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    ax=plt.gca()
    im=ax.imshow(np.transpose(zreion_table),interpolation='bilinear',
               origin='lower',extent=[np.log10(Rrange[0]),np.log10(Rrange[1]),deltaRrange[0],deltaRrange[1]])
    plt.xlabel(r'$\log\ R$ [Mpc]',fontsize=20)
    plt.ylabel(r'$\delta_{R}$',fontsize=20)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb=plt.colorbar(im,cax=cax)
    cb.set_label(r'$z_{\rm reion}$',fontsize=20)
    _cont = ax.contour(np.transpose(zreion_table),levels=np.linspace(4,20, 10), colors='white', alpha=0.5,
               origin='lower',extent=[np.log10(Rrange[0]),np.log10(Rrange[1]),deltaRrange[0],deltaRrange[1]])
    cb.add_lines(_cont)
    #plt.show()
    plt.savefig('zreion2.png',bbox_inches='tight')
    print('zreion should increase towards higher R and delta_R')
    print('this is not an exhaustive test, just a single quick check')

#I don't think we neeed the plots in here???
#added them again for test

#Look at zreion table

#plot_zreion_from_table(zreion_tablefile)

#show_zreion_from_table(zreion_tablefile)


