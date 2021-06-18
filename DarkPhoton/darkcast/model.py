# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import os, sys, inspect, math, collections
from . import pars

###############################################################################
# Update the model paths.
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "models"))
if os.getenv("DARKCAST_MODEL_PATH"):
    for path in reversed(os.getenv("DARKCAST_MODEL_PATH").split(":")):
        sys.path.insert(1, path)

###############################################################################
class ModelError(Exception):
    """
    Simple exception for the 'Model' class.
    """
    pass

###############################################################################
class Model:
    """
    Provides the information and methods needed to define a given
    model, e.g. 'dark_photon'.

    name:   name of the model.
    xfs:    dictionary of fermion couplings. Each coupling is a function 
            dependent upon mass (GeV).
    q:      quark U(3) charge matrix.
    """
    ###########################################################################
    def __init__(self, name, states = None, iwidth = None, path = None):
        """
        Load a model, given its name.

        The model must exist in the form '<name>.py' and is searched
        for along these paths in the following order:
        (0) The current directory within the Python interpreter.
        (1) The paths defined by the environment variable 
            'DARKCAST_MODEL_PATH'.
        (2) The 'models' directory of the Darkcast package.

        Each model must contain a fermion coupling dictionary named
        'xfs', where each coupling can either be a constant, or a mass 
        dependent function.

        The list 'states' may be defined, specifying the allowed final
        states for the model, e.g. ['e_e', 'mu_mu', 'invisible',
        ...]. Only these final states are used when calculating the
        total width. If not defined, all visible and invisible final
        states are used when calculating the total width.

        Optionally, an 'iwidth' function provides the invisible width
        for the model, given a mass and model and taking the form
        'iwidth(mass (GeV), model)'. Consequently, the invisible width
        can be defined as a function of the visible width. If no
        'iwidth' is defined, the invisible width is taken as zero. The
        invisible width is assumed to be dependent on the square of
        the global coupling. See the example model for further
        details.
        
        name:   name of the model.
        states: optionally, specify the allowed final states of the model.
        iwidth: optionally, specify the invisible width as a function of 
                a given mass and this model.
        path:   optionally, specify the path to load the module from.
        """
        # Import the model.
        self.name = name
        self.__cache = {}
        if path: sys.path.insert(1,path)
        model = __import__(name)
        if path: del sys.path[1]

        # Load the model's fermion couplings.
        self.xfs = {}
        for f in pars.mfs:
            try:
                float(model.xfs[f])
                self.xfs[f] = lambda m, f = f: float(model.xfs[f])
            except: 
                try: self.xfs[f] = model.xfs[f]
                except: raise ModelError(
                        "Error loading '%s' coupling from '%s'." % (f, name))

        # Load the model's invisible width function.
        try: self.__iwidth = iwidth if iwidth != None else model.iwidth
        except: self.__iwidth = lambda m, model: 0.0
        self.__iwidth(0, self)

        # Create the quark U(3) charge matrix.
        self.q = [self.xfs["u"], self.xfs["d"], self.xfs["s"]]

        # Load the model's defined final states.
        try: self.__states = states if states != None else model.states
        except: self.__states = ["visible", "invisible"]
        self.width("total", 0)
        try: self.width("total", 0)
        except: raise ModelError(
            "Invalid definition of allowed final states from '%s'." % name)

    ###########################################################################
    def trq(self, m, t):
        """
        Return the trace of the quark U(3)-charge matrix for the model
        with the diagonal of a given matrix, e.g. a meson generator T.
        
        m: mass at which to evaulate the couplings (GeV).
        t: diagonal of the matrix to perform the trace with, must be
           size 3.
        """
        try: return (t[0]*self.xfs["u"](m) + t[1]*self.xfs["d"](m) +
                     t[2]*self.xfs["s"](m))
        except: raise ModelError(
            "Invalid diagonal provided to the trace.")

    ###########################################################################
    def width(self, states, m, g = 1.0):
        """
        Return the width, in GeV, for the specified states, mass,
        and global coupling.

        states: final state or states.
        m:      mass (GeV).
        g:      global coupling (unitless).
        """
        # Loop over the states.
        total = 0
        for state in (states,) if isinstance(states, str) else states:

            # Use cached result if valid.
            cache = self.__cache.get(state)
            if cache and cache[0] == m: total += cache[-1]; continue
    
            # Invisible, visible, hadronic, and total widths.
            dtrs = state.split("_")
            if state == "invisible":
                part = self.__iwidth(m, self)
            elif state == "visible":
                part = self.width(
                    ["e_e", "mu_mu", "tau_tau", "nue_nue", "numu_numu", 
                     "nutau_nutau", "c_c", "b_b", "t_t", "hadrons"], m)
            elif state == "hadrons":
                part = self.width(pars.rfs.keys(), m)
            elif state == "total":
                part = self.width(self.__states, m)
    
            # Perturbative decay into a fermion pair, equation 2.13.
            elif len(dtrs) == 2 and dtrs[0] == dtrs[1] and dtrs[0] in pars.mfs:
                dtr = dtrs[0]
                cf, mf, xf = pars.cfs[dtr], pars.mfs[dtr], self.xfs[dtr](m)
                if m > 2.0*mf: part = (cf*xf**2.0*m/(12.0*math.pi)*(
                        1.0 + 2.0*mf**2/m**2)*math.sqrt(1.0 - 4.0*mf**2.0/m**2))
                else: part = 0

            # Decay into hadrons, equations 2.17 and 2.18.
            elif state in pars.rfs:
                part = 0
                for mesons, rf in pars.rfs[state].items():
                    sub = 1
                    for meson in mesons:
                        sub *= pars.rvs[meson]*self.trq(m, pars.tms[meson])
                    sub *= sub if len(mesons) == 1 else 2
                    sub *= rf(m)
                    part += m/(12*math.pi)*sub

            else: raise ModelError(
                "Unknown state '%s'." % state)

            # Cache the result.
            total += part
            self.__cache[state] = (m, part)
        return g*g*total

    ###########################################################################
    def tau(self, m, g = 1.0):
        """
        Return the lifetime, in seconds, for the specified mass and
        and global coupling.

        m: mass (GeV).
        g: global coupling (unitless).
        """
        return pars.hbar/self.width("total", m, g)

    ###########################################################################
    def g(self, m, tau):
        """
        Return the global coupling, for the specified mass and lifetime.

        m:   mass (GeV).
        tau: lifetime (seconds).
        """
        return math.sqrt(self.tau(m)/tau)

    ###########################################################################
    def bfrac(self, states, m):
        """
        Return the branching fraction for the specified states and mass.

        states: final state or states.
        m:      mass (GeV).
        """
        num = self.width(states, m)
        if num == 0: return 0.0
        den = self.width("total", m)
        if den == 0: return 0.0
        return num/den

###############################################################################
class Models(collections.OrderedDict):
    """
    Loads all 'Model's along the provided paths. The 'Models' object
    acts as an ordered dictionary for the individual models.
    """
    ###########################################################################
    def __init__(self, paths = None, states = None, iwidth = None):
        """
        Load all available models along the specified paths.

        paths:  paths to search for models. If no paths are specified,
                search the paths specified by DARKCAST_MODEL_PATH and
                the local Darkcast model directory.
        states: optionally, specify the allowed final states of the models.
        iwidth: optionally, specify the invisible width as a function of 
                a given mass and model.
        """
        super(Models, self).__init__()

        # Set default search paths.
        if paths == None:
            paths = []
            if os.getenv("DARKCAST_MODEL_PATH"): paths += [
                p for p in os.getenv("DARKCAST_MODEL_PATH").split(":")]
            paths += [os.path.join(os.path.dirname(os.path.realpath(
                            inspect.getfile(inspect.currentframe()))),"models")]

        # Load the models.
        for path in (paths,) if not hasattr(paths, "__iter__") else paths:
            models = sorted(os.listdir(path))
            for model in models:
                if not model.endswith(".py"): continue
                try: self[model[0:-3]] = Model(
                        model[0:-3], states, iwidth, path)
                except: pass
