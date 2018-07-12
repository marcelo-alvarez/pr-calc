import math
import numpy as np
import sys

import radialprofile

def powerspectrum(f):

    ff = np.fft.fftshift(np.fft.fft2(f))

    ff = ff / ff.shape[0] / ff.shape[1]

    ps2 = np.abs(ff)**2

    r, ps = radialprofile.azimuthalAverage(ps2, returnradii = True)

    return np.column_stack((r, ps))    

