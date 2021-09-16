#!/usr/bin/env python
# coding: utf-8

# Imports

# In[1]:


import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import peaks
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.integrate import quad
from colossus.lss import mass_function
import matplotlib.pyplot as plt


# Variables + cosmology

# In[2]:


#we set our own cosmology
#params         = {'flat': True, 'H0': 69.0, 'Om0': 0.2863, 'Ob0': 0.0463, 'sigma8': 0.82, 'ns': 0.96}  #Check this with pr-calc
params         = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.96}
cosmo          = cosmology.setCosmology('myCosmo', params)

h              = cosmo.H0 / 100.

# Mmin and zeta
Mmin           = 3E9 / h # in Msun (M in Colossus is in Msun/h) --------- This doesn't make M the same as in pr-calc anymore
Mmax           = 1E15/ h  
zeta           = 100  #70

# create grid of delta_R and R and the zreion table
# can mess with the grid later but an initial guess follows
nR             = 32
ndelta_R       = 33 # odd ndelta_R for a delta_R=0 entry with symmetric range
Rvalsmin       = 1.0
Rvalsmax       = 1e3
delta_Rmin     = -1
delta_Rmax     =  1
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
Rmin_Mpc       = (3*Mmin/4/np.pi/rho_Mpc)**(1./3.)   #R=(3M/4 pi)**1/3
Rmin           = Rmin_Mpc*h # Rmin is used by colussus and is in Mpc/h (multiply by )

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

# In[3]:


# function for fcoll = PS conditional collapsed fraction (ccfPS)
def ccfPS(z, R, delta_R):
    delta_c   = peaks.collapseOverdensity(corrections = True,z = z)
    sigma_R   = cosmo.sigma(R, z)  #R in Mpc/h
    sigma_min = cosmo.sigma(Rmin, z) #Rmin in Mpc/h
    arg       = ((delta_c - delta_R) / (np.sqrt(2.*sigma_min**2 - sigma_R**2)))
    f         = erfc(arg)
    return f
ccfPS = np.vectorize(ccfPS)


# Function to compute unconditional fcoll from colossus for a given hmf model and mass definition

# In[4]:


def integrand(lnM, rho_Mpc, redshift, model, mdef):
    '''
    This function computes the collapse fraction as an integral of the mass function indicated by the user
    '''
    h  = cosmo.H0 / 100. 
    M  = np.exp(lnM)     # input mass assumed to be in Msun
    M *= h               # convert from Msun to Msun/h for Colussus
    
    dndlnM = mass_function.massFunction(M, redshift, mdef = mdef, model = model, q_out = 'dndlnM')

    dndlnM *= h**3 # convert mass function from 1/(Mpc/h)^3 to 1/Msun/Mpc^3
    M      /= h    # convert mass back to Msun from Msun/h 
    
    dfcolldlnM = 1/rho_Mpc * M * dndlnM   # convert dfcoll/dlnM = 1/rho*M*dn/dlnM

    return dfcolldlnM

def ucf(z,hmf_model,hmf_mdef):
    # integrate hmf to obtain collapsed fraction 
    from scipy import integrate
    uncond_fcoll,   err = quad(integrand, np.log(Mmin), np.log(Mmax), 
                                args=(rho_Mpc, z, hmf_model, hmf_mdef))
    return uncond_fcoll
ucf = np.vectorize(ucf)


# Function to compute conditional fcoll using equation 3 from overleaf

# In[5]:


def ccf(z, R, delta_R, ucf_func, ucfPS_func):
    # eqn 3 from overleaf 
    # cond_fcoll(z,R,delta_R) = uncond_fcoll(z) / uncond_fcollPS * cond_fcollPS
    return ucf_func(z) / ucfPS_func(z) * ccfPS(z, R, delta_R) 


# Function to compute reionization redshift at a given radius and density contrast

# In[6]:


def zreion(R, delta_R, ucf_func, ucfPS_func):
    # note zvals is a global variable!!!
    fcoll_zvals = ccf(zvals,R,delta_R,ucf_func, ucfPS_func) 
    zeta_times_fcoll_zvals = zeta * fcoll_zvals
    z_of_zeta_times_fcoll = interp1d(zeta_times_fcoll_zvals, zvals, fill_value="extrapolate")
    zreion = z_of_zeta_times_fcoll(1.0)
    return zreion
zreion = np.vectorize(zreion)


# Populate collapsed fraction tables and set up collapsed fraction interpolation tables

# In[7]:


# tables for fcoll(z)
#   these are the unconditional collapsed fractions for PS 
#   as well as that for the specified hmf model
uncond_fcoll   = ucf(zvals,hmf_model,hmf_mdef).astype(np.float32)
uncond_fcollPS = ucf(zvals,'press74','fof').astype(np.float32)
ucf_interp     = interp1d(zvals, uncond_fcoll)
ucfPS_interp   = interp1d(zvals, uncond_fcollPS)


# Populate zreion table

# In[8]:


# table for zreion(R,delta_R)
#   this is the reionization redshift at a given overdensity and scale 
#   from the condition
#   zeta * fcoll(R,delta_R,zreion) = 1
zreion_table = zreion(Rvals2d, delta_Rvals2d, ucf_interp, ucfPS_interp).astype(np.float32)



#for i in range(len(zetavals)):
#    zreion_table = zreion(Rvals2d, delta_Rvals2d, ucf_interp, ucfPS_interp, zetavals[i]).astype(np.float32)


# Write tables in binaries with minimal headers easily readable in C

# In[9]:


# fcoll(z) (unconditional mass function)  ......why writing the unconditional to file?
fcoll_tablefile='fcoll_Mmin'+str(Mmin)+str(hmf_model)+'.tab'
f = open(fcoll_tablefile,'wb')
zrange       = np.asarray([zvals[0],zvals[-1]]).astype(np.float32)
dim          = np.asarray([len(zvals)]).astype(np.int32)
fcoll_table  = uncond_fcoll.astype(np.float32)
zrange.tofile(f)
dim.tofile(f)
fcoll_table.tofile(f)
f.close()

# zreion(R,delta_R) (reionization redshifts) 
#for i in range(len(zetavals)):
zreion_tablefile='zreion_zeta'+str(zeta)+'_Mmin'+str(Mmin)+'_'+str(hmf_model)+'.tab'
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

# In[10]:


# function to read zreion table
def read_zreiontable(fname):
    f = open(fname,'rb')
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
    #plt.savefig('zreion.png',bbox_inches='tight')
    plt.show()
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
    plt.show()
    print('zreion should increase towards higher R and delta_R')
    print('this is not an exhaustive test, just a single quick check')


# Look at zreion table

# In[11]:


#for i in range(len(zetavals)):
plot_zreion_from_table(zreion_tablefile)


# In[12]:


#for i in range(len(zetavals)):
show_zreion_from_table(zreion_tablefile)


# Check that uncond_fcoll == cond_fcoll when r is large and dR ~0

# In[13]:


plt.plot(zvals, ucfPS_interp(zvals), 
         'k-', label = 'Colussus integration of press74') #
plt.plot(zvals, ccfPS(zvals,Rvals[-1],delta_Rvals[dRzero_index]), 
         'r:', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=20)
plt.savefig('conditional.png',bbox_inches='tight')
plt.show()

print('red dotted line should match black line')
print('this is not an exhaustive test, just a single quick check')


# In[14]:


#samee plot as above but multiplied by zeta
plt.plot(zvals, zeta*ucfPS_interp(zvals), 
         'k-', label = 'Colussus integration of press74') #
plt.plot(zvals, zeta*ccfPS(zvals,Rvals[-1],delta_Rvals[dRzero_index]), 
         'r:', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=20)
plt.savefig('conditional.png',bbox_inches='tight')
plt.show()

print('red dotted line should match black line')
print('this is not an exhaustive test, just a single quick check')


# In[15]:


# print(len(zreion_table))
# mass_vec = np.logspace(np.log10(Mmin), np.log10(Mmax), 32)
# print(len(mass_vec))
# plt.semilogy(zreion_table, mass_vec)
# plt.show()


# In[16]:


#plot_zreion_from_table(test.zreion)
# file2 = open("test.zreion", "rb")
# bdata. decode(‘ASCII’)
# print(file2.read(20))
# file2.close()


# In[17]:


# Rrange       = np.fromfile(file2,dtype=np.float32, count=2)
# deltaRrange  = np.fromfile(file2,dtype=np.float32, count=2)
# dims         =     np.fromfile(file2,dtype=np.int32,   count=2)
# zreion = np.fromfile(file2,dtype=np.float32, count=dims[0]*dims[1])
# zreion = np.reshape(file2,(dims[0],dims[1]))


# In[18]:


#plots of fcoll

plt.plot(zvals, ccfPS(zvals,Rvals[0],delta_Rvals[dRzero_index]), 
         'r', label='conditional EPS for R='+str(Rvals[0])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]
plt.plot(zvals, ccfPS(zvals,Rvals[10],delta_Rvals[dRzero_index]), 
         'b', label='conditional EPS for R='+str(Rvals[10])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]
plt.plot(zvals, ccfPS(zvals,Rvals[-1],delta_Rvals[dRzero_index]), 
         'g', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.show()


# In[19]:


#plots of fcoll

plt.plot(zvals, ccfPS(zvals,Rvals[-1],delta_Rvals[0]), 
         'r', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[0])) # zvals, rvals,deltarvals   #[:,9,27]
plt.plot(zvals, ccfPS(zvals,Rvals[-1],delta_Rvals[dRzero_index]), 
         'b', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[dRzero_index])) # zvals, rvals,deltarvals   #[:,9,27]
plt.plot(zvals, ccfPS(zvals,Rvals[-1],delta_Rvals[-1]), 
         'g', label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[-1])) # zvals, rvals,deltarvals   #[:,9,27]

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=20)
plt.show()


# In[20]:


#plots of fcoll
for i in range(len(delta_Rvals)):
    plt.plot(zvals, ccfPS(zvals,Rvals[-1], delta_Rvals[i])
             , label='conditional EPS for R='+str(Rvals[-1])+', and deltaR='+str(delta_Rvals[i]))
    #plt.plot(zvals, ccfPS(zvals,Rvals[0], delta_Rvals[i]), linestyle=':'
             #, label='conditional EPS for R='+str(Rvals[0])+', and deltaR='+str(delta_Rvals[i]))


plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=10)
plt.ylim(1e-5,1)
plt.show()


# In[21]:


#plots of fcoll
for i in range(len(Rvals)):
    plt.plot(zvals, ccfPS(zvals, Rvals[i], delta_Rvals[-1])
             , label='conditional EPS for R='+str(Rvals[i])+', and deltaR='+str(delta_Rvals[-1]))

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.yscale('log')
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=8)
plt.show()


# In[22]:


# check table read and zreion
def plot_zreion_from_table(fname):
    #plot z_reion 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
    Rrange, deltaRrange, dims, zreion_table = read_zreiontable(fname)
    Rvals      = np.logspace(np.log10(Rrange[0]),     np.log10(Rrange[1]),     dims[0])
    deltaRvals = np.logspace(         deltaRrange[0],          deltaRrange[1], dims[1])
    plt.plot(delta_Rvals[:-1], zreion_table[:,0]           ,'r:',label='$R=$'+str(Rvals[0]))
    plt.plot(delta_Rvals[:-1], zreion_table[:,25]          ,'k-',label='$R=$'+str(Rvals[25]))
    plt.plot(delta_Rvals[:-1], zreion_table[:,-1]          ,'b:',label='$R=$'+str(Rvals[-1]))

    plt.legend(fontsize=20)
    plt.tick_params(labelsize=20)
    #plt.xscale('log')
    plt.xlabel('$\delta_R$',fontsize=20)
    plt.ylabel(r'$z_{\rm reion}$',fontsize=20)
    plt.savefig('zreion.png',bbox_inches='tight')
    plt.show()


# In[23]:


plot_zreion_from_table(zreion_tablefile)


# In[24]:


print(len(zreion_table))


# In[25]:


#plot_zreion_from_table(zreion_zeta100_Mmin4285714285.71_press74.tab)

f50         = open('zreion_zeta50_Mmin4285714285.71_press74.tab','rb')
Rrange50       = np.fromfile(f50,dtype=np.float32, count=2)
deltaRrange50  = np.fromfile(f50,dtype=np.float32, count=2)
dims50         = np.fromfile(f50,dtype=np.int32,   count=2)
zreion_table50 = np.fromfile(f50,dtype=np.float32, count=dims50[0]*dims50[1])
zreion_table50 = np.reshape(zreion_table50,(dims50[0],dims50[1]))

f75         = open('zreion_zeta75_Mmin4285714285.71_press74.tab','rb')
Rrange75       = np.fromfile(f75,dtype=np.float32, count=2)
deltaRrange75  = np.fromfile(f75,dtype=np.float32, count=2)
dims75         = np.fromfile(f75,dtype=np.int32,   count=2)
zreion_table75 = np.fromfile(f75,dtype=np.float32, count=dims75[0]*dims75[1])
zreion_table75 = np.reshape(zreion_table75,(dims75[0],dims75[1]))

f100         = open('zreion_zeta100_Mmin4285714285.71_press74.tab','rb')
Rrange100       = np.fromfile(f100,dtype=np.float32, count=2)
deltaRrange100  = np.fromfile(f100,dtype=np.float32, count=2)
dims100         = np.fromfile(f100,dtype=np.int32,   count=2)
zreion_table100 = np.fromfile(f100,dtype=np.float32, count=dims100[0]*dims100[1])
zreion_table100 = np.reshape(zreion_table100,(dims100[0],dims100[1]))


# In[26]:


import pandas as pd

pr_calc50 = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/HaloMassFunction/test50.history', sep='\s+', header=None)
pr_calc50.columns = ['z', 'chi', 'fcoll']

pr_calc75 = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/HaloMassFunction/test75.history', sep='\s+', header=None)
pr_calc75.columns = ['z', 'chi', 'fcoll']

pr_calc100 = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/HaloMassFunction/test100.history', sep='\s+', header=None)
pr_calc100.columns = ['z', 'chi', 'fcoll']

plt.plot(pr_calc50['z'], pr_calc50['chi'], 'r', label='pr_calc, zeta=50', linewidth =4)
plt.plot(pr_calc75['z'], pr_calc75['chi'], 'b', label='pr_calc, zeta = 75', linewidth =4)
plt.plot(pr_calc100['z'], pr_calc100['chi'], 'g', label='pr_calc, zeta = 100', linewidth =4)

plt.axvline(5.6, color='r', linestyle='-', linewidth = 1, label = 'ionization fraction = 1')
plt.axvline(6.1, color='b', linestyle='-', linewidth = 1)
plt.axvline(6.4, color='g', linestyle='-', linewidth = 1)

plt.axvline(zreion_table50[-1:,dRzero_index], linestyle='--', color='r', label='zreion from tables at R = 1000 and dR=0.0')
plt.axvline(zreion_table75[-1:,dRzero_index], linestyle='--', color='b')#, label='zeta=75')
plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')
plt.axhline(0.5, linestyle=':')
plt.axhline(0.9, linestyle=':')

#plt.axvline(zreion_table[0:,dRzero_index], linestyle=':')
plt.xlim(5,12)
plt.xlabel('redshift, z')
plt.ylabel('Ionization fraction')
plt.legend()
plt.show()

#print(zreion_table[:,-1])
#print(zreion_table[:,dRzero_index])
#print(delta_Rvals[dRzero_index])


# In[27]:


print(pr_calc100[87180:87205])


# In[28]:


pr_calc100 = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/HaloMassFunction/test100.history', sep='\s+', header=None)
pr_calc100.columns = ['z', 'chi', 'fcoll']

plt.semilogy(pr_calc100['z'], pr_calc100['fcoll'], 'g', label='pr_calc, zeta = 100', linewidth =4)

#plt.axvline(6.399066, color='g', linestyle='-', linewidth = 1)

#plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')

#plt.axhline(0.5, linestyle=':')
#plt.axhline(0.9, linestyle=':')

plt.xlim(0,20)
plt.xlabel('redshift, z')
plt.ylabel('Ionization fraction')
plt.legend()
plt.show()


# In[29]:




#plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')
#plt.axvline(6.4, color='g', linestyle='-', linewidth = 1)

plt.plot(pr_calc50['z'], pr_calc50['chi'], 'r', label='pr_calc, zeta=50', linewidth =4)
plt.plot(pr_calc75['z'], pr_calc75['chi'], 'b', label='pr_calc, zeta = 75', linewidth =4)
plt.plot(pr_calc100['z'], pr_calc100['chi'], 'g', label='pr_calc, zeta = 100', linewidth =4)

fz50 = 50 * (pr_calc50['fcoll'])
fz75 = 75 * (pr_calc75['fcoll'])
fz100 = 100 * (pr_calc100['fcoll'])

plt.plot(pr_calc50['z'], fz50, 'r', label='zeta x fcoll')
plt.plot(pr_calc75['z'], fz75, 'b')#, label='zeta x fcoll')
plt.plot(pr_calc100['z'], fz100, 'g')#, label='zeta x fcoll')#/check1
#plt.axvline(check1[dRzero_index], color='k', label='fff')

plt.axvline(zreion_table50[-1:,dRzero_index], linestyle='--', color='r', label='zreion from tables at R = 1000 and dR=0.0')
plt.axvline(zreion_table75[-1:,dRzero_index], linestyle='--', color='b')#, label='zeta=75')
plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')

plt.xlim(5,12)
plt.ylim(-0.1,1.2)
plt.legend()
plt.show()

#print(len(check1))
#print(zeta)
#print(pr_calc50['fcoll'][:-5])
#print(check1[:-5])


# In[30]:


#unconditional fcoll from colossus
#unfz50 = 50 * ucfPS_interp(zvals)
#unfz75 = 75 * (pr_calc75['fcoll'])
unfz100 = 100 * ucfPS_interp(zvals)

plt.plot(zvals, unfz100, 'k', label='zeta x uncondfcoll', linewidth =8)

#plt.plot(pr_calc50['z'], pr_calc50['chi'], 'r', label='pr_calc, zeta=50', linewidth =4)
#plt.plot(pr_calc75['z'], pr_calc75['chi'], 'b', label='pr_calc, zeta = 75', linewidth =4)
plt.plot(pr_calc100['z'], pr_calc100['chi'], 'g', label='pr_calc, zeta = 100', linewidth =4)

#fcoll from pr-calc
fz50 = 50 * (pr_calc50['fcoll'])
fz75 = 75 * (pr_calc75['fcoll'])
fz100 = 100 * (pr_calc100['fcoll'])

#plt.plot(pr_calc50['z'], fz50, 'r', label='zeta x fcoll')
#plt.plot(pr_calc75['z'], fz75, 'b')
plt.plot(pr_calc100['z'], fz100, 'g')

#plt.axvline(zreion_table50[-1:,dRzero_index], linestyle='--', color='r', label='zreion from tables at R = 1000 and dR=0.0')
#plt.axvline(zreion_table75[-1:,dRzero_index], linestyle='--', color='b')#, label='zeta=75')
plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')

plt.axhline(1.0, linestyle=':')

plt.xlim(5,12)
plt.ylim(-0.1,1.2)
#plt.ylim(0.1,10)
#plt.yscale('log')
plt.legend()
plt.show()


# In[31]:


plt.semilogy(pr_calc100['z'], pr_calc100['fcoll'], 'g', label='pr_calc, zeta = 100')
plt.semilogy(zvals, ucfPS_interp(zvals), 'k', label='interpolation')

plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.axhline(1, linestyle=':',c='k')
plt.legend(fontsize=10)
plt.ylim(1e-5,1)
plt.xlim(0,15)
plt.show()


# In[ ]:


pr_calc100bigL = pd.read_csv('testlambda500zeta100.history', sep='\s+', header=None)
pr_calc100bigL.columns = ['z', 'chi', 'fcoll']

pr_calc100bigL2 = pd.read_csv('testlambda800zeta100.history', sep='\s+', header=None)
pr_calc100bigL2.columns = ['z', 'chi', 'fcoll']

pr_calc100bigL3 = pd.read_csv('testlambda1000zeta100.history', sep='\s+', header=None)
pr_calc100bigL3.columns = ['z', 'chi', 'fcoll']

pr_calc100bigL4 = pd.read_csv('testlambda1500zeta100.history', sep='\s+', header=None)
pr_calc100bigL4.columns = ['z', 'chi', 'fcoll']

plt.plot(pr_calc100['z'], pr_calc100['chi'], 'g', label='pr calc, zeta = 100, lambda = 300')
plt.plot(pr_calc100bigL['z'], pr_calc100bigL['chi'], 'k', label='pr calc, zeta = 100, lambda = 500')
plt.plot(pr_calc100bigL2['z'], pr_calc100bigL2['chi'], 'b', label='pr calc, zeta = 100, lambda = 800')
plt.plot(pr_calc100bigL3['z'], pr_calc100bigL3['chi'], 'y', label='pr calc, zeta = 100, lambda = 1000')
plt.plot(pr_calc100bigL4['z'], pr_calc100bigL4['chi'], 'r', label='pr calc, zeta = 100, lambda = 1500')

plt.axvline(zreion_table100[-1:,dRzero_index], linestyle='--', color='g')#, label='zeta=100')


plt.xlim(5,12)
plt.ylim(-0.1,1.2)
plt.legend()
plt.show()


# In[ ]:


print(len(unfz50))
print(len(pr_calc50['z']))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def show_zreion_from_table(fname):
    #plot z_reion 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
    Rrange, deltaRrange, dims, zreion_table = read_zreiontable(fname)
    Rvals      = np.logspace(np.log10(Rrange[0]),     np.log10(Rrange[1]),     dims[0])
    deltaRvals = np.logspace(     deltaRrange[0],          deltaRrange[1],     dims[1])
    
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
    plt.show()


# In[ ]:


show_zreion_from_table(zreion_tablefile)


# In[ ]:


# check table read and zreion
def plot_zreion_from_table2(fname):
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
    plt.show()
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
    plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


import pandas as pd

pr_calc = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/cmbonly_spectra_dr4.01/test.history', sep='\s+', header=None)
pr_calc.columns = ['z', 'chi', 'fcoll']

plt.plot(pr_calc['z'], pr_calc['chi'], 'r', label='pr_calc')
plt.axvline(6.0)
for i in range(len(zetavals)):
    plt.axvline(zreion_table[-1:,dRzero_index], linestyle=':')

#plt.axvline(zreion_table[0:,dRzero_index], linestyle=':')
plt.xlim(5,12)
plt.xlabel('redshift, z')
plt.ylabel('Ionization fraction')
plt.show()

#print(zreion_table[:,-1])
print(zreion_table[:,dRzero_index])
print(delta_Rvals[dRzero_index])



# In[ ]:





# In[ ]:





# In[ ]:


print(delta_Rvals)


# In[ ]:


print(len(zvals))


# In[ ]:


#checking zreion with pr-cals
import pandas as pd

zreion = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/HaloMassFunction/test.history', sep='\s+', header=None)
zreion.columns = ['z', 'chi', 'fcoll', 'zreion']

plt.plot(zreion['z'], zreion['chi'], label='pr_calc')


# In[ ]:




