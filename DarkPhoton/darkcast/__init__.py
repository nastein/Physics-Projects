# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
"""
Darkcast is the companion software package to the paper 'Serendipity
in dark photon searches' and is a framework for recasting constraints
from dark photon searches into other models. The following directory
structure is in place:

vmd:    contains the R_mu^f interpolation grids from equation 2.5, split 
        by individual meson contributions, including interference.
models: all models reported in the paper are defined here.
limits: limits reported in the paper are provided here. Files of the 
        form '<name>.py' define all needed information for a limit,
        e.g. meta-data, model, production, bounds, and
        efficiency. Files of the form '<name>.lmt' provide the actual
        limit data, while' <name>.prd' provide production mechanism
        ratios for more complex scenarios.

Example usage is as follows:
'''
# Load the module.
import darkcast

# Change any global parameters, here the speed of light (m/s).
darkcast.pars.c = 3e8

# Load a limit, here the LHCb prompt limit.
limit = darkcast.Limit('LHCb_Aaij2017rft_prompt')

# Print the notes and BibTex for the limit.
print limit.notes
print limit.bibtex

# Load a model for recasting.
model = darkcast.Model('B_boson')

# Recast from the limit model to the new model.
recast = limit.recast(model)

# Write out the recast limit.
recast.write('LHCb_B_boson.txt')
'''

More detailed examples with explanations are provided in the
'examples' directory and further documentation is provided per
sub-module and class.
"""
from .utils import Dataset, Datasets
from .model import Model, Models
from .production import BreitWigner, Production
from .efficiency import Efficiency
from .limit import Limit, Limits
