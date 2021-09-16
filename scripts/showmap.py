#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import array
import math
import numpy
import matplotlib.pyplot as plt
import sys

#taufile = open('z_4096_t0.06_l032_m1e8.tmap',"rb")
taufile = open(sys.argv[1],"rb")

n=array.array('i')
fov=array.array('f')
tau=array.array('f')

n.fromfile(taufile,2)
fov.fromfile(taufile,2)

print('dimensions of map are:  ',n[0],'x',n[1])
print('field of view is:       ',fov[0]/2./math.pi*360.,'x',fov[1]/2./math.pi*360.,' degrees')

tau.fromfile(taufile,n[0]*n[1])

taubar=sum(tau)/n[0]/n[1]

tauimg = numpy.zeros((n[0],n[1]))

for i in range(0,n[0]):
    for j in range(0,n[1]):
        k = i*n[1] + j
        tauimg[i][j]=tau[k]
        tau[k]=tau[k]**2    

tau=numpy.reshape(tau,(n[0]*n[1]))

print('mean of tau is:         ',taubar)
print('variance of tau is:     ',sum(tau)/n[0]/n[1]-taubar**2)
print('rms in percent:         ',100.*math.sqrt(sum(tau)/n[0]/n[1]-taubar**2)/taubar)

v = numpy.linspace(0.08,0.1,11, endpoint=True)
imgplot = plt.imshow(tauimg)
imgplot.set_clim(0.08,0.1)
plt.colorbar(ticks=v)
#plt.show()
plt.savefig(sys.argv[2],format='pdf')
