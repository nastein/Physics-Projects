# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# This is an example user limit. To include new limits within
# Darkcast, please contact the authors or create a merge request at
# https://gitlab.com/philten/darkcast/.

# Import the Darkcast module. Note that the path does not need to
# specified. Since the Darkcast module must already have been imported
# to load this limit the module must exist on the path.
import darkcast

###############################################################################
# Notes on the limit can be given with the 'notes' variable. This does
# not have to be provided, and a default note is set.
notes = """
This is an example limit which demonstrates how new limits can be
written.
"""

###############################################################################
# BibTeX information for the limit can also be provided, taken
# directly from http://inspirehep.net/. Again, BibTeX information is
# optional.
bibtex = """
@article{Ilten:2018crw,
 author         = "Ilten, Philip and Soreq, Yotam and Williams, Mike and
                   Xue, Wei",
 title          = "{Serendipity in dark photon searches}",
 year           = "2018",
 eprint         = "1801.04847",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "MIT-CTP/4976, CERN-TH-2017-282",
 SLACcitation   = "%%CITATION = ARXIV:1801.04847;%%"
}
"""

###############################################################################
# A model for the limit must be defined and must be of the type
# 'darkcast.Model'. All the provided limits are defined with the dark
# photon model but could be defined using other models.
model = darkcast.Model("dark_photon")

###############################################################################
# The new boson production mechanism(s) must also be defined, and be
# of the type 'darkcast.Production'. The simplest production
# definition is a single mechanism which can be from the following:
#
#        p_brem: proton-beam bremsstrahlung.
#        A_brem: A-beam bremsstrahlung, where A can be any fundamental fermion.
#        A_A:    Drell-Yan, where A can be any fundamental fermion.
#        A:      vector meson mixing, where A is the vector meson.
#        A_B:    meson decays of the form A -> B + X, where A is the decaying 
#                meson, B is the SM daughter and X is the NP daughter.
#
# The following defines production from eta decays, eta -> gamma X:
#
# production = darkcast.Production("eta_gamma")

# A custom mechanism may also be passed, where the mechanism must take
# the X boson mass and a model as arguments. The mechanism is assumed
# to be dependent upon the square of the global coupling. Here,
# proton-beam bremsstrahlung is defined (note this does not use the mass).
#
# def p_brem(m, model): return (2*model.xfs['u'] + model.xfs['d'])**2
# production = darkcast.Production(p_brem)

# It is also possible to define production from multiple mechanisms by
# passing a dictionary where the keys are the mechanisms (either
# pre-defined or user) and the values are the relative fraction of
# that mechanism. See CHARM_Bergsma1985qz for an example. Here,
# production is defined as 1% eta -> gamma X and 99% omega -> pi0 X:
#
# production = darkcast.Production({"eta_gamma": 0.01, "omega_pi0": 0.99})

# The mechanism fractions can be mass dependent functions. Note that for
# each mass the fractions are not normalized. It is the responsibility of
# the user to ensure normalized fractions.
#
# def eta_frac(m): return m**2/(m**2 + 5.0)
# def omega_frac(m): return 5.0/(m**2 + 5.0)
# production = darkcast.Production({"eta_gamma": eta_frac,
#                                   "omega_pi0": omega_frac})

# Data grids using the 'Datasets' class may be used to define the
# production dictionary, see for example LHCb_Aaij2017rft_prompt. The
# first column of the dataset must be the X boson mass (in GeV) and
# all remaining columns each specify the fraction for a given
# mechanism. The first line of the file must specify the production
# mechanism for each column, separated by spaces. For example the
# first line:
#
# # mass eta_gamma omega_pi0
#
# would specify the second column as eta -> gamma X production and the
# third column as omega -> pi0 X production.  Linear interpolation is
# used between mass points, and frozen at the edge of the dataset.
#
# Datasets are searched for along the following paths:
# (0) The absolute path, if the absolute path is given.
# (1) The current directory within the Python interpreter.
# (2) The paths defined by the environment variable 'DARKCAST_DATA_PATH'.
# (3) The Darkcast package directory.
# Below, the dataset is being found with path (1).
production = darkcast.Production(darkcast.Datasets("user_limit.prd"))

###############################################################################
# Next, the decay for the limit must be specified. This can be any
# valid final state(s) available for the 'Model' class. The following
# sets a di-pion decay. See BaBar_Lees2014xha for an example with
# multiple final states.
decay = "pi+_pi-"

###############################################################################
# The bounds for the limit must also be defined, through a data grid
# loaded with either the 'Dataset' or 'Datasets' classes. If the bound
# is a lower bound, then the 'Dataset' class is used. Only two columns
# should be provided in the data file, the first column is the X boson
# mass in GeV and the second column is the bound on the coupling.
#
# bounds = darkcast.Dataset("user_limit_single.lmt")
#
# An example of this type of bound is given by APEX_Abrahamyan2011gv.

# If the full bound is provided, e.g. r-values as is done for
# LHCb_Aaij2017rft_displaced, a 'Datasets' object should also be
# used. Here, the first column is still mass in GeV, the second column
# is the global coupling, and the third column is the r-value. Each
# mass point must have the same global coupling points defined.
#
# bounds = darkcast.Dataset("user_limit_rvalue.lmt")

# If a lower and upper bound is defined, e.g. E137_Bjorken1988as, then
# the 'Datasets' class is used. The data grid must have three columns:
# the X boson mass in GeV, the lower bound on the global coupling, and
# the upper bound on the global coupling. The first line of the file
# must read:
#
# # mass lower upper
#
# The bound is then defined as follows:
#
# bounds = darkcast.Datasets("user_limit_double.lmt")

###############################################################################
# Finally, an efficiency needs to be defined. This is done by
# specifying a proper time fiducial with the values t0 and t1 and
# defining the efficiency as e^(-t0/tau) - e^(-t1/tau).

# For a prompt one-sided limit t0 will typically be 0, and if the
# efficiency is always unity, than t1 is infinity.
#
# efficiency = darkcast.Efficiency(t0 = 0, t1 = float('inf'))

# For full limits with bounds defined by r-values, the efficiency
# should be initialized as:
#
# efficiency = darkcast.Efficiency(rval = True)
#
# Outside the valid bound range, the upper efficiency is defined by t0
# = tau_min and t1 = infinity, and the lower efficiency is defined by
# t0 = 0 and t1 = tau_max. This type of efficiency may only be used
# with full r-value limits.

# If the limit is two-sided, the values for t0 and t1 can determined
# by solving:
#
# g_min^2 efficiency(tau_min) = g_max^2 efficiency(tau_max)
#
# where t1 = t0(1 + L_decay/L_shield) such that L_decay is the length
# of the decay volume and L_shield is the length of the shielding
# volume. Only the length ratio needs to be provided. This type of
# efficiency may only be used with a double-sided limit.
#
# efficiency = darkcast.Efficiency(lratio = 100.0/1.0)

###############################################################################
# For convenience, the above bounds/efficiency combinations are
# collected below.

# Single-sided lower limit.
# bounds = darkcast.Dataset("user_limit_single.lmt")
# efficiency = darkcast.Efficiency(t0 = 0, t1 = float('inf'))

# Double-sided lower/upper limit.
# bounds = darkcast.Datasets("user_limit_double.lmt")
# efficiency = darkcast.Efficiency(lratio = 100.0/1.0)

# Full r-value limit.
bounds = darkcast.Dataset("user_limit_rvalue.lmt")
efficiency = darkcast.Efficiency(rvals = True)
