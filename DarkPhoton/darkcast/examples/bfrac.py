# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

# This example calculates the branching fractions for every model
# available in darkcast/models. If matplotlib is available these
# branching fractions are then plotted. The produced figures
# correspond to figure 3 of Ilten:2018crw.

# Update the system path to find the Darkcast module.
# This assumes that 'examples' is in 'darkcast/examples.'
import sys, os, inspect, itertools
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "../../"))

# Import the Darkcast module.
import darkcast

# Load all the available models in darkcast/models and DARKCAST_MODEL_PATH.
models = darkcast.Models()

# Alternatively, all the models from a folder, '/foo/bar', could be
# loaded as:
#
# models = darkcast.Models("/foo/bar")
#
# Note that any models in the directory that are not valid will not be
# loaded. A single model with name 'foo_bar' can be loaded as:
#
# model = darkcast.Model("foo_bar")

# Create the dictionary of channels for which to calculate the
# branching fractions. A channel can be either a single state,
# e.g. 'mu_mu' or a list of final states, e.g. ['mu_mu', 'e_e']. The
# following can be commented in or out depending on the channels
# needed. Note this is an 'OrderedDict' to ensure that the keys remain
# in the order specified. Keys which require mathmode in LaTeX are
# enclosed in '$'.
import collections
channels = collections.OrderedDict([
        # Entries take the form (key, value).
        
        # All available fundamental fermion pairs.
        ("$e_e$",                 "e_e"),
        ("$mu_mu$",               "mu_mu"),
        #("$tau_tau$",             "tau_tau"),
        #("$nue_nue$",             "nue_nue"),
        #("$numu_numu$",           "numu_numu"),
        #("$nutau_nutau$",         "nutau_nutau"),
        #("$d_d$",                 "d_d"), # Included in exclusive hadrons.
        #("$u_u$",                 "u_u"), # Included in exclusive hadrons.
        #("$s_s$",                 "s_s"), # Included in exclusive hadrons.
        #("$c_c$",                 "c_c"),
        #("$b_b$",                 "b_b"),
        #("$t_t$",                 "t_t"),

        # Combine Neutrinos into a single channel.
        ("$nu_nu$",               ["nue_nue", "numu_numu", "nutau_nutau"]),

        # All available exclusive hadronic states.
        #("$pi+_pi-$",             "pi+_pi-"),
        #("$pi+_pi-_pi+_pi-$",     "pi+_pi-_pi+_pi-"),
        #("$pi+_pi-_pi0_pi0$",     "pi+_pi-_pi0_pi0"),
        #("$pi+_pi-_pi0$",         "pi+_pi-_pi0"),
        #("$pi0_gamma$",           "pi0_gamma"),
        #("$K_K$",                 "K_K"),
        #("$K_K_pi$",              "K_K_pi"),
        #("other hadrons",         "other"),

        # Alias for all exclusive hadronic final states above.
        ("hadrons",                "hadrons"),

        # Alias for all visible final states, e.g. everything above
        # except 'd_d', 'u_u', and 's_s'.
        #("visible",                "visible"),

        # All invisible final states.
        #("invisible",              "invisible"),
        
        # All possible final states to consider when calculating the
        # total width. Typically 'visible' and 'invisible' but this
        # can be specified by the user when creating the model with
        # the 'states' variable, e.g. darkcast.Model('dark_photon',
        # states = ['e_e', 'mu_mu']).
        #("total",                  "total"),
        ])

# Create the list of masses for which to calculate the branching fractions.
masses = [mass*1e-2 for mass in range(1, 200)]

# Try to load matplotlib.
try: import matplotlib.pyplot as pyplot
except: pyplot = None
colors = ["red", "green", "blue", "orange", "magenta", "cyan", "gray"]

# Loop over all the models.
for name, model in models.items():

    # If possible, initialize the plot.
    if pyplot:
        fig, ax = pyplot.subplots()
        icolor, labels = itertools.cycle(colors), {}

    # Loop over the channels.
    for label, channel in channels.items():
        
        # Calculate the branching fraction for the model and channel
        # as a function of mass.
        bfracs = [model.bfrac(channel, mass) for mass in masses]

        # Additionally, the width can be calculated using the 'width'
        # method and the same channels as for 'bfrac'.
        #
        # widths = [model.width(channel, mass, g = 1) for mass in masses]

        # Save the branching fraction to a text file.
        txt = open("bfrac_%s_%s.txt" % (name, label.replace("$", "")), "w")
        for mass, bfrac in zip(masses, bfracs):
            txt.write("%11.4e %11.4e\n" % (mass, bfrac))
        txt.close()
            
        # Plot. The 'latex' utility converts commonly used symbols to
        # LaTeX, e.g. 'pi0' -> '\pi^{0}'.
        if pyplot:
            if not label in labels: color = next(icolor); labels[label] = color
            else: color = labels[label]; label = None
            ax.plot(masses, bfracs, label = darkcast.utils.latex(label),
                    color = color)

    # Save the plot.
    if pyplot:
        legend = ax.legend(loc = "best", fontsize = 10)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 1])
        ax.set_xlabel("mass [GeV]")
        ax.set_ylabel("branching fraction")
        ax.set_title(darkcast.utils.latex(model.name))
        darkcast.utils.logo()
        fig.savefig("bfrac_%s.pdf" % name)
