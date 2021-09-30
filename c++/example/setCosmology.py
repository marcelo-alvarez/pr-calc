#!/usr/bin/env python

#we set our own cosmology
from colossus.cosmology import cosmology
#def param_names():
#    params         = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.96}
#    cosmo          = cosmology.setCosmology('myCosmo', params)
#    h              = cosmo.H0 / 100.
#    return cosmo

params         = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.96}
cosmo          = cosmology.setCosmology('myCosmo', params)
h              = cosmo.H0 / 100.
print(h)

