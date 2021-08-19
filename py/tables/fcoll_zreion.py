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
params        = {'flat': True, 'H0': 69.0, 'Om0': 0.282, 'Ob0': 0.042, 'sigma8': 0.81, 'ns': 0.96}  #Check this with pr-calc
cosmo         = cosmology.setCosmology('myCosmo', params)

h             = cosmo.H0 / 100.

# Mmin and zeta
Mmin          = 3E9 / h # in Msun
Mmax          = 1E15/ h  #1E12/ h
zeta          = 50

# create grid of delta_R and R and the zreion table
# can mess with the grid later but an initial guess follows
nR            = 10
ndelta_R      = 50
Rvalsmin      = 1.0
Rvalsmax      = 1e3
delta_Rmin    = -2
delta_Rmax    =  2
Rvals         = np.logspace(np.log10(Rvalsmin), np.log10(Rvalsmax), nR)
delta_Rvals   = np.linspace(delta_Rmin, delta_Rmax, ndelta_R)
zreion_table  = np.zeros((nR, ndelta_R))

# create grid of z and the fcoll table
dz            = 0.1
zmin          = 0.0 # just for completeness
zmax          = 20  # most reionization models we don't care about z>20 
nz            = int((zmax - zmin) / dz)
fcoll_table   = np.zeros(nz)

rho_Mpc       = 2.775E11 * cosmo.Om0 * h**2
 
zvals         = np.linspace(zmin, zmax, nz) # zvals will be approximately spaced at dz
fcoll_table   = np.zeros(nz)

# convert Mmin to Rmin
Rmin          = (3*Mmin/4/np.pi/rho_Mpc)**(1./3.)

#print(Rvals)
#print(zvals)
print(delta_Rvals)


# In[3]:


#print(zvals)
#print(delta_Rvals)
#print(Rvals)

print(len(zvals))
print(len(delta_Rvals))
print(len(Rvals))


# Function to compute PS conditional fcoll

# In[4]:


# function for fcoll conditional
def cond_fcoll(z, delta_R, R):
    delta_c   = peaks.collapseOverdensity(corrections=True,z=z)
    #sigma_min = cosmo.sigma(Rmin,z)
    sigma_R   = cosmo.sigma(R,z)
    arg       = ((delta_c - delta_R) / (np.sqrt(2.*sigma_min**2 - sigma_R**2)))
    f         = erfc(arg)
    return f

cond_fcollPS = np.zeros((len(zvals), len(Rvals), len(delta_Rvals)))
for i in range(len(zvals)):
    #print(i,nz)
    sigma_min    = cosmo.sigma(Rmin,zvals[i])
    delta_c      = peaks.collapseOverdensity(corrections = True, z = zvals[i])
    for j in range(len(Rvals)):
        sigma_R      = cosmo.sigma(Rvals[j], zvals[i], j=0, filt='sharp-k', inverse=False, derivative=False, kmin=None, kmax=None, ps_args={'model': 'eisenstein98', 'path': None})
        for k in range(len(delta_Rvals)): 
            cond_fcollPS[i,j,k] = cond_fcoll(zvals[i], delta_Rvals[k], Rvals[j])


# In[5]:


#print(cond_fcollPS)
#print(cond_fcollPS[:,:,:])#same as above, all of fcoll corresponding to all of zvals, rvals and delta_Rvals
#print(cond_fcollPS[0,:,:])# prints fccoll at z=0, for all of rvals and deltarvals. zvalls index goes from 0 t0 199 (because it's len =200) which corresponds to z 20
#print(cond_fcollPS[:,0,:])#prints out all the fccoll at the first value of rvals (in this casee 1), for all values of zvals and deltarvals. index goes to 10 which correspnds to rvals  1000
print(cond_fcollPS[:,:,0])#prints out all the fcoll at the first value of deltarvals (in this case -2), for all values of zvals and rvals. index goes to 10 which corresponds to deltarvals of 2
print(cond_fcollPS[:,:,0].shape)


# In[6]:


plt.plot(zvals, cond_fcollPS[:,9,9], 'r')
#plt.plot(cond_fcollPS[20,:,9], cond_fcollPS[0,0,:], 'r--')
plt.plot(zvals, cond_fcollPS[:,-5,5], 'b')
plt.plot(zvals, cond_fcollPS[:,4,7], '--g')
plt.yscale('log')

plt.xlabel('z', fontsize=20)
plt.ylabel('F_coll', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=10)
plt.show()

print(cond_fcollPS[20,:,9])  #Rvals


# making different plots to further analayse the data
# 

# In[7]:



plt.plot(zvals, cond_fcollPS[:,0,0], '--g')
plt.plot(zvals, cond_fcollPS[:,0,5], 'r')

plt.yscale('log')
plt.show()


# In[8]:


#plot fcoll at several different values of z, 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
plt.plot(cond_fcollPS[:,0,0])
plt.show()


# compute unconditional fcoll from colossus. Ideally, this would be done for PS nad ST

# In[43]:


def intergrand_fcoll(M, rho_Mpc, redshift, model):
    
    '''
    This function computes the collapse fraction as an inteegral of the mass function indicated by the user
    '''
    
    for model in mass_function.models:
        if 'fof' in mass_function.models[model].mdefs:
            mfunc  = mass_function.massFunction(M, redshift, model = model, mdef = 'fof', q_out = 'dndlnM')
        else:
            mfunc = mass_function.massFunction(M, redshift, mdef = '200m', model = model)
            
            M     /= 0.69
            mfunc  = mfunc * (1/M)   #convert dndlnM to dndM
            
            mfunc *= 0.69**3
            return (1/rho_Mpc) * mfunc * M 
        if M > (1e13 / 0.69):
            return 0
    

F_coll = np.zeros((len(zvals), len(Rvals), len(delta_Rvals)))
F_coll_PS = np.zeros((len(zvals), len(Rvals), len(delta_Rvals)))

#print(F_coll)

for i in range(len(zvals)):
    F_coll[i,0,0], err = quad(intergrand_fcoll, Mmin, Mmax, args=(rho_Mpc, zvals[i], 'press74'), epsabs = 1e-6, limit=100)
    F_coll_PS[i,0,0], err = quad(intergrand_fcoll, Mmin, Mmax, args=(rho_Mpc, zvals[i], 'press74'), epsabs = 1e-6, limit=100)
    


# In[44]:


print(F_coll[:,0,0])
#print(F_coll_PS)
print(1e-20)


# check that uunconditional fcoll aat a given reedshift is equal to conditionaal fcoll when R is large and dR ~ 0

# In[45]:


#check that uncond_fcoll aat a given z == cond_fcoll when r is large and dR ~0
plt.plot(zvals, cond_fcollPS[:,9,27], 'r:', label='conditional at R=1000, and deltaR ~0') # zvals, rvals,deltarvals
#plt.plot(zvals, cond_fcollPS[:,1,13], 'r', label='conditional at R=1000, and deltaR ~0') # zvals, rvals,deltarvals
plt.plot(zvals, F_coll_PS[:,0,0], '-', label = 'unconditional') #

plt.xlabel('z', fontsize=20)
plt.ylabel('F_coll', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=10)
plt.yscale('log')
plt.axhline(1, linestyle='--')
plt.legend()
plt.show()

#red dotted line should match blue curve


# compute the final conditional fcoll from equation 3 in overleeaf

# In[42]:



print(F_coll.shape)
print(cond_fcollPS.shape)

#F_coll_final = 

F_coll = np.nan_to_num(F_coll)
F_coll_PS = np.nan_to_num(F_coll_PS)
cond_fcollPS[:,:,:] = np.nan_to_num(cond_fcollPS[:,:,:])

F_coll_final = (F_coll * cond_fcollPS[:,:,:]) / F_coll_PS #eqn 3...F_coll_PS[i] = bar(fcollPS), = fcollPS
print(F_coll_final)

#print(F_coll_final.shape)
#print(F_coll_PS[1,:,:])


# plot the values to be sure they agree since we are begining with just PS for checks

# In[18]:


plt.plot(zvals, (F_coll_PS[:,0,0]), 'b-', label = 'unconditional from colosus') #


plt.semilogy(zvals,F_coll_final[:,0,0], 'r-', label='Final conditional fcoll from equation')
plt.xlabel('z', fontsize=20)
plt.ylabel('F_coll', fontsize=20)
    
plt.plot(zvals,cond_fcollPS[:,-1,5], 'g:', label='conditional fcoll PS')
plt.yscale('log')
plt.legend()
plt.show()

#red curve should be equal to greeen


# In[ ]:


# populate the tables
fcoll_table = F_coll_final
#print(fcoll_table)


# In[ ]:


for i in range(len(Rvals)):
    for j in range(len(delta_Rvals)):
        cond_fcoll_zvals = cond_fcoll(zvals, delta_Rvals[j], Rvals[i])
        zeta_times_fcoll_zvals = zeta * cond_fcoll_zvals
        z_of_zeta_times_fcoll = interp1d(zeta_times_fcoll_zvals, zvals, fill_value="extrapolate")
        zreion_table[i,j] = z_of_zeta_times_fcoll(1.0)
zreion_table = zreion_table.astype(np.float32)


# In[ ]:


print(len(zvals))
print(len(F_coll_final))
print(len(Rvals))


# In[ ]:


zreion_table = np.nan_to_num(zreion_table)
print(zreion_table[:,49])#all values of zreion at rvals index 2
print(zreion_table.shape)


# In[ ]:


#plot z_reion 1) vs R at a few values of delta_R, and 2) vs delta_R at a few values of R
plt.plot(Rvals, zreion_table[:,49], 'o')
plt.yscale('log')
plt.show( )


# In[ ]:


# write the tables in binaries with minimal headers easily readable
# in C
# fcoll table
#wrote tables in .txt as well so I can see the output

tablefile2='fcoll_Mmin'+str(Mmin)+'.txt'  #add the mass function model to the file name
#tablefile2='fcoll_Mmin'+str(Mmin)+str(mass_function.models[model])+'.txt'
f2 = open(tablefile2,'w')
fcoll_table.tofile(f2, sep=" ", format="%f")
np.savetxt(tablefile2, np.column_stack([zvals, fcoll_table[:,0,0]]), fmt='%.18g')#zvals, Rvals, delta_Rvals, fcoll_table

tablefile='fcoll_Mmin'+str(Mmin)+'.tab'
f = open(tablefile,'wb')
zrange = [zmin, zmax]
zrange = np.asarray(zrange).astype(np.float32)
dim   = [nz]
dim   = np.asarray(dim).astype(np.int32)
zrange.tofile(f)
dim.tofile(f)
fcoll_table.tofile(f)


# In[ ]:


# zreion table
tablefile3='zreion_zeta'+str(zeta)+'_Mmin'+str(Mmin)+'.txt'
#print(tablefile)
f3 = open(tablefile3,'w')
zreion_table.tofile(f3, sep=" ", format="%f")
np.savetxt(tablefile3, np.column_stack([Rvals, zreion_table]), fmt='%.18g')

tablefile='zreion_zeta'+str(zeta)+'_Mmin'+str(Mmin)+'.tab'
f = open(tablefile,'wb')
ranges = [Rvalsmin,Rvalsmax,delta_Rmin,delta_Rmax]
ranges = np.asarray(ranges).astype(np.float32)
dims   = [nR, ndelta_R]
dims   = np.asarray(dims).astype(np.int32)
ranges.tofile(f)
dims.tofile(f)
zreion_table.tofile(f)


# In[ ]:


f.close()
f2.close()
f3.close()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def intergrand_fcoll(M, rho_Mpc, redshift):
    mfunc  = mass_function.massFunction(M, redshift, mdef = 'fof', model = 'press74', q_out = 'dndlnM')
    mfunc  = mfunc * (1/M)   #convert dndlnM to dndM
    M     /= 0.69
    mfunc *= 0.69**3
    return (1/rho_Mpc) * mfunc * M    

F_coll = np.zeros(len(zvals))
F_coll_PS = np.zeros(len(zvals))

for i in range(len(zvals)):
    F_coll[i], err = quad(intergrand_fcoll, Mmin, Mmax, args=(rho_Mpc, zvals[i]))
    F_coll_PS[i], err = quad(intergrand_fcoll, Mmin, Mmax, args=(rho_Mpc, zvals[i]))
    

