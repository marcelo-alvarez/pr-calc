#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import array
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 6:
    print('usage: ./map2pk.py <input map> <fileroot> <min> <max> <field>\n')
    print(sys.argv)
    sys.exit(2)

mapfile = open(sys.argv[1],"rb")
outfilepdf = sys.argv[2]+".pdf"
outfilepng = sys.argv[2]+".png"
cmin    = float(sys.argv[3])
cmax    = float(sys.argv[4])
field   = sys.argv[5]


n=array.array('i')
fov=array.array('f')
mapdata=array.array('f')

# Read size and Field of view
n.fromfile(mapfile,2)
fov.fromfile(mapfile,2)

# Set limit to be half the size of map in degrees
limit = fov[0]/2/math.pi*360/2

# Read images, need to convert to float64 (double) to avoid precision issues (affects mean and var)

mapdata.fromfile(mapfile,n[0]*n[1])
img=np.zeros((n[0],n[1]))
for i in range(0,n[0]):
    for j in range(0,n[1]):
        k = i*n[1] + j
        img[i][j]=mapdata[k]

print('min, max of image: ',np.amin(img),np.amax(img))

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
print(outfilepdf)
plt.savefig(outfilepdf,format='pdf')
plt.savefig(outfilepng,format='png')

