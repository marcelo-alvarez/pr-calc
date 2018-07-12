#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 6:
    print '\nusage: ./map2pk.py <input map> <pdffile> <min> <max> <field>\n'
    print sys.argv
    sys.exit(2)

mapfile = open(sys.argv[1],"rb")
outfilepdf = open(sys.argv[2]+".pdf","w")
outfilepng = open(sys.argv[2]+".png","w")
cmin    = float(sys.argv[3])
cmax    = float(sys.argv[4])
field   = sys.argv[5]

# Read size and Field of view
n    = np.fromfile(mapfile, dtype=np.int32,   count=2)
fov = np.fromfile(mapfile, dtype=np.float32, count=2)

# Set limit to be half the size of map in degrees
limit = fov[0]/2/math.pi*360/2

# Read images, need to convert to float64 (double) to avoid precision issues (affects mean and var)

img = np.fromfile(mapfile, dtype=np.float32, count=n[0]*n[1]).astype(np.float64).reshape(n)
img = img[:n[0],:n[1]]

print np.amin(img),np.amax(img)

if field=='dtb':
    img/=1e3
    cmax/=1e3

imgplot = plt.imshow(img,interpolation='none')

imgplot.set_clim(cmin,cmax)
imgplot.set_extent([-limit,limit,-limit,limit])

plt.xlabel('x [degrees]')
plt.ylabel('y [degrees]')

#imgplot.axes.get_xaxis().set_visible(False)
#imgplot.axes.get_yaxis().set_visible(False)

cb = plt.colorbar()

if field=="ksz": cb.set_label(r'${\Delta}T_{\rm kSZ}$ [ ${\mu}$K ]')
if field=="dtb": cb.set_label(r'${\Delta}T_{\rm b}$ [ mK ]')
if field=="tau": cb.set_label(r'${\tau}_{\rm es}$')

plt.savefig(outfilepdf,format='pdf')
plt.savefig(outfilepng,format='png')

