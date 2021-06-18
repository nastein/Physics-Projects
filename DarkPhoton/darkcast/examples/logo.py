# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# This example just draws the Darkcast logo, in case anyone needs it!

# Update the system path to find the Darkcast module.
# This assumes that 'examples' is in 'darkcast/examples.'
import sys, os, inspect, itertools
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "../../"))

# Load the Darkcast module.
import darkcast

# Draw the logo.
from matplotlib import pyplot
darkcast.utils.logo(0, 0, 1)
pyplot.savefig("logo.pdf", bbox_inches = "tight")
