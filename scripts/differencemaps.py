#!/usr/bin/env python

import math
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print '\nusage: ./map2pk.py <input map 1> <input map 2> <output map>\n'
    print sys.argv
    sys.exit(2)

file1 = open(sys.argv[1],"rb")
file2 = open(sys.argv[2],"rb")
file3 = open(sys.argv[3],"wb")

# Read size and Field of view
n1   = np.fromfile(file1, dtype=np.int32,   count=2)
fov1 = np.fromfile(file1, dtype=np.float32, count=2)

n2   = np.fromfile(file2, dtype=np.int32,   count=2)
fov2 = np.fromfile(file2, dtype=np.float32, count=2)

# Read images, need to convert to float64 (double) to avoid precision issues (affects mean and var)

img1 = np.fromfile(file1, dtype=np.float32, count=n1[0]*n1[1])

img2 = np.fromfile(file2, dtype=np.float32, count=n2[0]*n2[1])

img3 = img1 - img2

n1.tofile(file3)
fov1.tofile(file3)
img3.tofile(file3)

