# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# This is an example user model. To include new models within
# Darkcast, please contact the authors or create a merge request at
# https://gitlab.com/philten/darkcast/.

###############################################################################
# The constants from Darkcast can be accessed as follows. Here,
# 'help(darkcast.pars)' will return the documentation on the relevant
# parameters. Note that changing any of the parameters here will
# change them globally.
#
# import darkcast          # Load Darkcast.
# help(darkcast.pars)      # Print the help on the parameters.
# print darkcast.pars.hbar # Print the value of hbar.

###############################################################################
# The fermion couplings must be defined via the 'xfs'
# dictionary. These couplings can be either a constant, or a function
# which takes a single argument of mass in GeV.
xfs = {
    "e":       0.0,
    "mu":      0.0,
    "tau":     0.0,
    "nue":     0.0,
    "numu":    0.0,
    "nutau":   0.0,
    "d":       lambda m: m**2,
    "u":      -1.0,
    "s":      20.0,
    "c":     -20.0,
    "b":       1.0,
    "t":      -1.0
    }

###############################################################################
# An invisible width can be optionally defined as 'iwidth'. This must
# be a function that takes as an argument a mass (in GeV) and the
# model itself. This means that the invisible width can be defined as
# a function of the visible width, which is done here. See
# BaBar_Lees2017lec as a limit where this is done.
def iwidth(m, model): return model.width("visible", m)

# The invisible width can also be set as just a number:
#
# iwidth = 0
# 
# By default, the invisible width is set as 0.

###############################################################################
# The final states to consider when calculating the total width, can
# also be optionally set. By default, all possible final states are
# considered (the exclusive hadronic final states, not the
# perturbative light quark components, are used). The 'states'
# variable can be either a list of final states, or just a single
# final state. Here the states are defined as the default: visible and
# invisible.
states = ["visible", "invisible"]
