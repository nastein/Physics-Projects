# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# This example demonstrates how both a new model and a new limit can
# be created. The example model is 'user_model.py' and the example
# limit is 'user_limit.py' with the associated files
# 'user_limit_single.lmt', 'user_limit_double.lmt',
# 'user_limit_rvals.lmt', and 'user_limit.prd'. To add a model or
# limit to the official Darkcast repository please contact the authors
# or create a merge request at https://gitlab.com/philten/darkcast/.

# Update the system path to find the Darkcast module.
# This assumes that 'examples' is in 'darkcast/examples.'
import sys, os, inspect, itertools
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "../../"))

# Import the Darkcast module.
import darkcast

# New models and limits must exist in the form '<name>.py' and are
# searched for along these paths in the following order:
# (0) The current directory within the Python interpreter.
# (1) The paths defined by the environment variables 
#     'DARKCAST_MODEL_PATH' and 'DARKCAST_LIMIT_PATH' respectively.
# (2) The 'models' and 'limits' directories of the Darkcast package, 
#     respectively.
# Since the example models exist in the current directory, they are
# picked up by path (0). The 'dark_photon' model is picked up by path
# (2).

# Load the new limit. 
userlimit = darkcast.Limit("user_limit")

# Load the new model.
usermodel = darkcast.Model("user_model")

# Try to load matplotlib.
try: import matplotlib.pyplot as pyplot
except: pyplot = None
colors = ["red", "green", "blue", "orange", "magenta", "cyan", "gray"]

# If possible, initialize the plot.
if pyplot:
    fig, ax = pyplot.subplots()
    icolor, labels = itertools.cycle(colors), {}

# Create the directory for recasted limits.
if not os.path.exists("recast"): os.makedirs("recast")

# Recast the user limit to the user model.
for model in [darkcast.Model("dark_photon"), usermodel]:

    # Recast the limit, this returns an object of type 'Datasets'.
    recast = userlimit.recast(model)
 
    # Save the limit to a text file. This is done with the
    # 'Datasets.write' method.
    if not os.path.exists("recast/" + model.name):
        os.makedirs("recast/" + model.name)
    recast.write("recast/%s/%s.lmt" % (model.name, userlimit.name))
        
    # Plot. The 'Datasets.plots' method returns formatted lists of
    # x and y points which can be easily passed to a plotting
    # package.
    if pyplot:
        label = darkcast.utils.latex(model.name)
        if not label in labels: color = next(icolor); labels[label] = color
        else: color = labels[label]; label = None
        for x, y in recast.plots():
            ax.fill(x, y, label = label, alpha = 0.3, color = color)
            label = None

# Save the plot.
if pyplot:
    legend = ax.legend(loc = "best", fontsize = 10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([2e-3, 1e2])
    ax.set_ylim([1e-8, 1e1])
    ax.set_xlabel("mass [GeV]")
    ax.set_ylabel("g")
    ax.set_title(darkcast.utils.latex(userlimit.name))
    darkcast.utils.logo()
    fig.savefig("user.pdf")
