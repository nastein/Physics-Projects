# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# Define the fermion couplings.
from darkcast.pars import ge
xfs = {
    "e":      -1*ge,
    "mu":     -1*ge,
    "tau":    -1*ge,
    "nue":     0*ge,
    "numu":    0*ge,
    "nutau":   0*ge,
    "d":      -1*ge/3,
    "u":       2*ge/3,
    "s":      -1*ge/3,
    "c":       2*ge/3,
    "b":      -1*ge/3,
    "t":       2*ge/3
    }
