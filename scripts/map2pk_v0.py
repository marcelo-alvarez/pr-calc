#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

import radialprofile
import powerspectrum

if len(sys.argv) != 2:
    print '\nusage: ./map2pk.py <input map> <output cl>\n'
    sys.exit(2)

taufile = open(sys.argv[1],"rb")
clfile  = open(sys.argv[2],"w")

# Read size and Field of view
n = np.fromfile(taufile, dtype=np.int32, count=2)
fov = np.fromfile(taufile, dtype=np.float32, count=2)

# Read image, need to convert to float64 (double) to avoid precision issues (affects mean and var)
tauimg = np.fromfile(taufile, dtype=np.float32, count=n[0]*n[1]).astype(np.float64).reshape(n)
tauimg = tauimg[:n[0],:n[1]]

print ''
print '***** Info for file:    ',sys.argv[1]
print 'dimensions of map are:  ',n[0],'x',n[1]
print 'field of view is:       ',fov[0]/2./math.pi*360.,'x',fov[1]/2./math.pi*360.,' degrees'

tau_mean = tauimg.mean()
tau_var  = tauimg.var()
tau_rms  = math.sqrt(tau_var)

print 'mean of tau is:         ', tau_mean
print 'variance of tau is:     ', tau_var
print 'rms of tau is:          ', tau_rms

scl = 2*math.pi / fov[0]

ti = tauimg - tauimg.mean()

tmps = powerspectrum(ti)

l  = tmps[:,0]
cl = tmps[:,1]

# convert l from map units [0:npixel] to multipoles [0:2pi/fov*npixel]
l = l * scl

# convert cl to 'spectral density'
cl = cl / (scl * scl) * (2*math.pi)*(2*math.pi)

# l*cl/2pi is power per l-mode
lcl = l*cl / (2*math.pi)

# l2cl is power per log-l, aka l^2*Cl/(2pi)
l2cl = l * lcl 

for i in range(0,tmps.shape[0]):
    print >>psfile, l[i], l2cl[i]

