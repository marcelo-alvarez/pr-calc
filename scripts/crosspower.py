import math
import numpy as np
import sys

import radialprofile

def crosspower(f1,f2):

    ff1 = np.fft.fftshift(np.fft.fft2(f1))
    ff2 = np.fft.fftshift(np.fft.fft2(f2))

    ff1 = ff1 / ff1.shape[0] / ff1.shape[1]
    ff2 = ff2 / ff2.shape[0] / ff2.shape[1]

    ps2 = np.real(ff1 * np.conjugate(ff2))

    r, ps = radialprofile.azimuthalAverage(ps2, returnradii = True)

    return np.column_stack((r, ps))    

