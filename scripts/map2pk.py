#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

import radialprofile
#import powerspectrum

from powerspectrum import powerspectrum

if len(sys.argv) != 3:
    print('\nusage: ./map2pk.py <input map> <output cl>\n')
    print(sys.argv)
    sys.exit(2)

taufile = open(sys.argv[1],"rb")
clfile  = open(sys.argv[2],"w")

# Read size and Field of view
n = np.fromfile(taufile, dtype=np.int32, count=2)
fov = np.fromfile(taufile, dtype=np.float32, count=2)

#n = np.array([512, 512], dtype="int32")
#fov = np.array([1.6 * np.pi / 180., 1.6 * np.pi / 180.], dtype="float32")

# Read image, need to convert to float64 (double) to avoid precision issues (affects mean and var)
tauimg = np.fromfile(taufile, dtype=np.float32, count=n[0]*n[1]).astype(np.float64).reshape(n)
tauimg = tauimg[:n[0],:n[1]]

print('')
print('***** Info for file:    ',sys.argv[1])
print('dimensions of map are:  ',n[0],'x',n[1])
print('field of view is:       ',fov[0]/2./math.pi*360.,'x',fov[1]/2./math.pi*360.,' degrees')

tau_mean = tauimg.mean()

ti = tauimg - tauimg.mean()

tau_var  = ti.var() 
tau_rms  = math.sqrt(tau_var)

print('mean of tau is:         ', tau_mean)
print('variance of tau is:     ', tau_var)

scl = 2*math.pi / fov[0]

tmps = powerspectrum(ti)

l  = tmps[:,0]
cl = tmps[:,1]

# convert l from map units [0:npixel] to multipoles [0:2pi/fov*npixel]
l = l * scl

cl = fov[0]**2 * cl

clll = (2 * l + 1) * cl / 4. / np.pi 

print("--->", np.trapz(clll,l), tau_var)

# l*cl/2pi is power per l-mode
lcl = l*cl / (2*math.pi)

# l^2*Cl/(2pi) is power per log-l
l2cl = l * lcl 

for i in range(0,tmps.shape[0]):
    #print >>clfile, l[i], clll[i], l2cl[i]
    print(l[i], clll[i], l2cl[i], file=clfile)

