# Cosmological parameters
Omegam        0.27
Omegal        0.73
Omegab        0.044
h             0.7
ns            0.96
sigma8        0.8
w             -1
	
# Box size, dimension, and redshift of input density field

BoxSize         4e3 # comoving box size in Mpc/h
N               512   #2048 # dimension of the grid
zInit           0
# This is the Reionization redshift file location.
# The file should be the reionization field

DeltaFile	./ICs/delta

# These are the parameters describing the map of tau

fov	    	20  # field of view in degrees
theta		30   # polar angle of center of map in degrees
phi		30    # azimuthal angle of center of map in degrees
NPixels		2048  # number of pixels per dimension in map
NSide		1024 # Healpix N_side parameter

# The parameters are for the actual raytracing

NRedshifts       10000
InitialRedshift  5.5
FinalRedshift    15
