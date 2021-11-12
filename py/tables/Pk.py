#!/usr/bin/env python

# coding: utf-8

# Imports
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import peaks
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.integrate import quad
from colossus.lss import mass_function
import matplotlib.pyplot as plt
import pickle
import sys

''' This file uses the current cosmology to compute the matter power spectrum and is run as python Pk.py param.col '''

with open(sys.argv[1], 'rb') as handle:
    b = pickle.load(handle)

# Variables + cosmology
#we set our own cosmology
params         = {'flat': b["flat"], 'H0': float(b["H0"]), 'Om0': float(b["Om0"]), 'Ob0': float(b["Ob0"]), 'sigma8': float(b["sigma8"]), 'ns': float(b["ns"])}
cosmo          = cosmology.setCosmology('myCosmo', params)
h              = cosmo.H0 / 100.

#Writing P(k) in tables
k = 10**np.arange(-5,2,0.02)  #wavenumber k (in comoving h/Mpc)
Pk = cosmo.matterPowerSpectrum(k, z=0.0, model='eisenstein98', path=None, derivative=False)
k = k * h
Pk = Pk / h**3

pkfile = np.column_stack([k, Pk])
np.savetxt("../../c++/example/ICs/pkfile.txt", pkfile, fmt=['%1.7e','%1.7e'])
