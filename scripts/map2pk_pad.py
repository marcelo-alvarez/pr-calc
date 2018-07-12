#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

import radialprofile
#import powerspectrum

from powerspectrum import powerspectrum

if len(sys.argv) != 3:
    print '\nusage: ./map2pk.py <input map> <output cl>\n'
    print sys.argv
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

# Find mean and variance
tau_mean = tauimg.mean()
tauimg = tauimg - tauimg.mean()
tau_var  = tauimg.var() 
tau_rms  = math.sqrt(tau_var)

print 'mean of tau is:         ', tau_mean
print 'variance of tau is:     ', tau_var

# Create window function
window = 0*tauimg+1

# Pad with zeros
fov=2*fov
zn=np.zeros(n[0])
z2n=np.zeros(2*n[0])

tauimg=np.column_stack((tauimg,[zn for i in range(n[0])]))
tauimg=np.row_stack((tauimg,[z2n for i in range(n[0])]))
tauimg=np.roll(tauimg,n[0],axis=0)
tauimg=np.roll(tauimg,n[0],axis=1)

window=np.column_stack((window,[zn for i in range(n[0])]))
window=np.row_stack((window,[z2n for i in range(n[0])]))
window=np.roll(window,n[0],axis=0)
window=np.roll(window,n[0],axis=1)

# Calculate power spectrum
tps = powerspectrum(tauimg)
wps = powerspectrum(window)

l  = tps[:,0]

cl  = tps[:,1]
clw = wps[:,1]

# divide by window function power spectrum
cl = cl / clw
#cl = cl * 4

# convert l from map units [0:npixel] to multipoles [0:2pi/fov*npixel]
scl = 2*math.pi / fov[0]
l = l * scl

cl = fov[0]**2 * cl

clll = (2 * l + 1) * cl / 4. / np.pi 

print "--->", np.trapz(clll,l), tau_var

# l*cl/2pi is power per l-mode
lcl = l*cl / (2*math.pi)

# l^2*Cl/(2pi) is power per log-l
l2cl = l * lcl 

for i in range(0,tps.shape[0]):
    print >>clfile, l[i], clll[i], l2cl[i]

