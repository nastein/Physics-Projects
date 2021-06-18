# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# Define the fermion couplings.
from darkcast.pars import ge
from math import pi
xfs = {
    "e":      -ge**2/(4*pi)**2,
    "mu":     -ge**2/(4*pi)**2,
    "tau":    -ge**2/(4*pi)**2,
    "nue":     0,
    "numu":    0,
    "nutau":   0,
    "d":       1./3,
    "u":       1./3,
    "s":       1./3,
    "c":       1./3,
    "b":       1./3,
    "t":       1./3
    }
