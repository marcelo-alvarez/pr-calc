#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

import radialprofile
#import powerspectrum

from powerspectrum import powerspectrum
from crosspower    import crosspower

if len(sys.argv) != 4:
    print '\nusage: ./map2pk.py <input map 1> <input map 2> <output cl>\n'
    print sys.argv
    sys.exit(2)

file1 = open(sys.argv[1],"rb")
file2 = open(sys.argv[2],"rb")
clfile = open(sys.argv[3],"w")

# Read size and Field of view
n1   = np.fromfile(file1, dtype=np.int32,   count=2)
fov1 = np.fromfile(file1, dtype=np.float32, count=2)

n2   = np.fromfile(file2, dtype=np.int32,   count=2)
fov2 = np.fromfile(file2, dtype=np.float32, count=2)

# Read images, need to convert to float64 (double) to avoid precision issues (affects mean and var)

img1 = np.fromfile(file1, dtype=np.float32, count=n1[0]*n1[1]).astype(np.float64).reshape(n1)
img1 = img1[:n1[0],:n1[1]]

img2 = np.fromfile(file2, dtype=np.float32, count=n2[0]*n2[1]).astype(np.float64).reshape(n2)
img2 = img2[:n2[0],:n2[1]]

mean1 = img1.mean()
mean2 = img2.mean()

img1 = img1 - mean1
img2 = img2 - mean2

var1  = img1.var() 
var2  = img2.var() 

rms1  = math.sqrt(var1)
rms2  = math.sqrt(var2)

print ''
print '***** Info for file:    ',sys.argv[1]
print 'dimensions of map are:  ',n1[0],'x',n1[1]
print 'field of view is:       ',fov1[0]/2./math.pi*360.,'x',fov1[1]/2./math.pi*360.,' degrees'
print 'mean:                   ',mean1
print 'variance:               ',var1

print ''
print '***** Info for file:    ',sys.argv[2]
print 'dimensions of map are:  ',n2[0],'x',n2[1]
print 'field of view is:       ',fov2[0]/2./math.pi*360.,'x',fov2[1]/2./math.pi*360.,' degrees'
print 'mean:                   ',mean2
print 'variance:               ',var2

scl1 = 2*math.pi / fov1[0]
scl2 = 2*math.pi / fov2[0]

ps1 = powerspectrum(img1)
ps2 = powerspectrum(img2)

cc  = crosspower(img1,img2)

l  = ps1[:,0]

cl1 = ps1[:,1]
cl2 = ps2[:,1]
cl  = cc[:,1]

# convert l from map units [0:npixel] to multipoles [0:2pi/fov*npixel]
l = l * scl1

cl1 = fov1[0]**2 * cl1
cl2 = fov1[0]**2 * cl2
cl  = fov1[0]**2 * cl

clll1 = (2 * l + 1) * cl1 / 4. / np.pi 
clll2 = (2 * l + 1) * cl2 / 4. / np.pi 
clll  = (2 * l + 1) * cl  / 4. / np.pi 

print "--->", np.trapz(clll1,l), var1
print "--->", np.trapz(clll2,l), var2

# l*cl/2pi is power per l-mode
lcl1 = l*cl1 / (2*math.pi)
lcl2 = l*cl2 / (2*math.pi)
lcl  = l*cl  / (2*math.pi)

# l^2*Cl/(2pi) is power per log-l
l2cl1 = l * lcl1 
l2cl2 = l * lcl2 
l2cl  = l * lcl 

for i in range(0,ps1.shape[0]):
    print >>clfile, l[i], l2cl1[i], l2cl2[i], l2cl[i]

