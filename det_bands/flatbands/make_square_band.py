"""
make_square_band.py

Usage:  nu_low nu_high delta_nu

writes two columns (nu, transmission) to standard output
(you can redirect to a file)

output nu runs from 0.9*nu_low to 1.1*nu_high

"""

import sys
import numpy as np

nu_low = float(sys.argv[1])
nu_high = float(sys.argv[2])
delta_nu = float(sys.argv[3])

freqs = np.arange(0.9*nu_low, 1.1*nu_high, delta_nu)
trans = np.zeros(len(freqs))
in_band = np.argwhere((freqs >= nu_low) & (freqs <= nu_high))
trans[in_band] = 1.0

A = np.stack((freqs,trans)).T

np.savetxt(sys.stdout.buffer, A, fmt='%.3f')

