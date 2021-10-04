#!/usr/bin/env python

# coding: utf-8

# Imports

#import sys
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
#params         = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.96}
params         = {'flat': b["flat"], 'H0': float(b["H0"]), 'Om0': float(b["Om0"]), 'Ob0': float(b["Ob0"]), 'sigma8': float(b["sigma8"]), 'ns': float(b["ns"])}
cosmo          = cosmology.setCosmology('myCosmo', params)
h              = cosmo.H0 / 100.

#Writing P(k) in tables
#matterPowerSpectrum(self, k, z=0.0, model='eisenstein98', path=None, derivative=False)

k = 10**np.arange(-5,2,0.02) / h  #wavenumber k (in comoving h/Mpc)
Pk = cosmo.matterPowerSpectrum(k, z=0.0, model='eisenstein98', path=None, derivative=False)
Pk *= h**3

pk_tablefile='pk.tab'
f = open(pk_tablefile,'wb')
k       = np.asarray(k).astype(np.float32)
Pk          = np.asarray(Pk).astype(np.int32)
k.tofile(f)
Pk.tofile(f)
f.close()
