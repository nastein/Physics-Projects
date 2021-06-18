# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import os, sys, inspect, collections
from . import utils

###############################################################################
# Update the limit paths.
sys.path.insert(1, os.path.join(os.path.dirname(os.path.realpath(
                inspect.getfile(inspect.currentframe()))), "limits"))
if os.getenv("DARKCAST_LIMIT_PATH"):
    for path in reversed(os.getenv("DARKCAST_LIMIT_PATH").split(":")):
        sys.path.insert(1, path)

###############################################################################
class LimitError(Exception):
    """
    Simple exception for the 'Limit' class.
    """
    pass

###############################################################################
class Limit:
    """
    Represents all the needed information to define a limit. This
    class contains the following members.
    
    notes:      text providing any relevant notes for the limit. This is
                optional, and if not defined a default will be assigned.
    bibtex:     the optional BibTex entry from inSPIRE for this limit.
    model:      the model defining the limit, must be of class 'Model'.
    production: production which contains the production mechanism or 
                mechanisms for the limit. Must be of class 'Production'.
    bounds:     the actual bounds for the limit. This can either be
                a 'Dataset' when a single-sided bound or 'Datasets' when
                a double-sided bound.
    efficiency: efficiency of class 'Efficiency'.
    """
    def __init__(self, name, path = None):
        """
        Initialize a limit, given its name.

        The limit must exist in the form <name>.py and is searched for
        along these paths in the following order:
        (0) The current directory within the Python interpreter.
        (1) The paths defined by the environment variable 
            'DARKCAST_LIMIT_PATH'.
        (2) The 'limits' directory of the Darkcast package.

        Each limit must have a model, production, decay, bounds, and
        efficiency defined. Optionally, notes and a BibTex entry can
        be provided.

        name: name of the limit.
        """
        # Check if the limit should be initialized.
        self.name = name

        # Import the limit.
        if path: sys.path.insert(1, path)
        limit = __import__(name)
        if path: del sys.path[1]

        # Load the notes and BibTeX.
        try: self.notes = limit.notes
        except: self.notes = "The limit %s has no notes." % name
        try: self.bibtex = limit.bibtex
        except: self.bibtex = "The limit %s has no BibTeX entry." % name

        # Load the model.
        if not hasattr(limit, "model"): raise LimitError(
            "No model is defined for '%s'." % name)
        self.model = limit.model
        self.model.width("total", 1)
        
        # Load the production.
        if not hasattr(limit, "production"): raise LimitError(
            "No production is defined for '%s'." % name)
        self.production = limit.production
        self.production.ratio(1, 1, 1, self.model, self.model)

        # Load the decay.
        try: self.decay = limit.decay
        except: raise LimitError(
            "No decay is defined for '%s'." % name)
        try: self.model.width(self.decay, 1)
        except: raise LimitError(
            "The decay defined for '%s' is not valid." % name)

        # Load the bounds.
        if not hasattr(limit, "bounds"): raise LimitError(
            "No bounds are defined for '%s'." % name)
        try:
            self.bounds = utils.Datasets()
            self.bounds["rvals" if limit.bounds.dim() == 2 
                        else "lower"] = limit.bounds
        except: self.bounds = limit.bounds
        for key, bound in self.bounds.items():
            if not isinstance(bound, utils.Dataset): raise LimitError(
                "The bound %s is not of type Dataset." % key)

        # Load the efficiency.
        if not hasattr(limit, "efficiency"): raise LimitError(
            "No efficiency is defined for '%s'." % name)
        self.efficiency = limit.efficiency
        self.efficiency.ratio(1, self, 1.0, 1.0)

    ###########################################################################
    def recast(self, model, gmax = 1e5):
        """
        Recast these limits to a given model. Returns a dictionary
        with entries of 'lower' and when relevant, 'upper'. Each entry
        is of the form: [[m0, m1, ...], [g0, g1, ...]].

        model: model for recasting, must be of type 'Model', e.g. 
               Model('dark_photon').
        gmax:  maximum coupling to recast. If a coupling is greater than 
               or equal to this, then that mass point is skipped.
        """
        # Initialize the recast bounds.
        rvals = self.bounds.get("rvals")
        lower = utils.Dataset()
        upper = utils.Dataset() if rvals or "upper" in self.bounds else None

        # Recast r-value bounds.
        if rvals: 

            # Loop over the r-values.
            jdx, tl, tu, gl, gu, rl, ru,  = 0, 0, 0, gmax, -gmax, 1, 1
            tmin, tmax = float("inf"), 0
            for idx, r1 in enumerate(rvals.vals):

                # Determine the r-value via equation 2.21.
                m, g1 = rvals._Dataset__f2x(idx)
                tau = self.model.tau(m, g1)
                g0 = model.g(m, tau)
                b0 = model.bfrac(self.decay, m)
                b1 = self.model.bfrac(self.decay, m)
                pr = self.production.ratio(m, g0, g1, model, self.model)
                r0 = r1*b1/(b0*pr) if b0*pr != 0 else gmax

                # Update the lower/upper limits.
                if r0 < 1:
                    if g0 < gl: tl, gl, rl = tau, g0, r0
                    if g0 > gu: tu, gu, ru = tau, g0, r0
                tmin = min(tmin, tau)
                tmax = max(tmax, tau)
                
                # Check if final g for given m.
                jdx += 1
                if jdx == len(rvals.axes[1]):

                    # Solve for any limits outside bounds.
                    if tl == tmax and rl < 1:
                        f = lambda g: g**2/gl**2*self.efficiency.ratio(
                            m, self, model.tau(m, g), tl)/rl - 1
                        try: gl = utils.solve(f, x1 = gl)
                        except: pass
                    if tu == tmin and ru < 1:
                        f = lambda g: g**2/gu**2*self.efficiency.ratio(
                            m, self, model.tau(m, g), tu)/ru - 1
                        try: gu = utils.solve(f, x0 = gu)
                        except: pass

                    # Update the bounds.
                    lower.axes[0].append(m); lower.vals.append(abs(gl))
                    upper.axes[0].append(m); upper.vals.append(abs(gu))
                    jdx, gl, gu, rl, ru = 0, gmax, -gmax, 1, 1
                    tmin, tmax = float("inf"), 0
            
        # Recast lower/upper bounds.
        else:

            # Loop over the masses.
            g0l = None
            for m, g1l in zip(self.bounds["lower"].axes[0],
                              self.bounds["lower"].vals):
    
                # Check if branching fraction and production is non-zero.
                if model.bfrac(self.decay, m) == 0: continue
                if self.model.bfrac(self.decay, m) == 0: continue
                if self.production.ratio(m, 1, 1, model, self.model) == 0:
                    continue
    
                # Check if limit is above maximum and set guess.
                if g1l >= gmax: 
                    lower.axes[0].append(m); lower.vals.append(gmax)
                    if upper: upper.axes[0].append(m); upper.vals.append(gmax)
                    continue
                ggl = g0l if g0l != None else g1l
    
                # Solve equation 2.2 for the lower bound.
                f = lambda g: (
                    model.bfrac(self.decay, m)/self.model.bfrac(self.decay, m)
                    * self.production.ratio(m, g, g1l, model, self.model) 
                    * self.efficiency.ratio(m, self, model.tau(m, g), 
                                            self.model.tau(m, g1l)) - 1)
                try: g0l = utils.solve(f, x = ggl)
                except: g0l = gmax
    
                # If upper bound, find second zero, e.g. equations C.5 - C.7.
                if upper == None:
                    lower.axes[0].append(m); lower.vals.append(g0l)
                else:
                    if g0l == gmax: g0u = gmax
                    elif f(g0l*1.01) < 0:
                        try: g0u = utils.solve(f, x1 = g0l*0.99)
                        except: g0u = gmax
                    else:
                        try: g0u = utils.solve(f, x0 = g0l*1.01)
                        except: g0u = gmax
                    if g0u < g0l: g0u, g0l = g0l, g0u
                    lower.axes[0].append(m); lower.vals.append(g0l)
                    upper.axes[0].append(m); upper.vals.append(g0u)
                                   
        # Return the recast bounds.
        bounds = utils.Datasets()
        bounds["lower"] = lower
        if upper != None: bounds["upper"] = upper
        return bounds

###############################################################################
class Limits(collections.OrderedDict):
    """
    Loads all 'Limit's along the provided paths. The 'Limits' object
    acts as an ordered dictionary for the individual models.
    """
    ###########################################################################
    def __init__(self, paths = None):
        """
        Load all available models along the specified paths.

        paths: paths to search for models. If no paths are specified,
               search the paths specified by DARKCAST_LIMIT_PATH and
               the local Darkcast limit directory.
        """
        super(Limits, self).__init__()

        # Set default search paths.
        if paths == None:
            paths = []
            if os.getenv("DARKCAST_LIMIT_PATH"): paths += [
                p for p in os.getenv("DARKCAST_LIMIT_PATH").split(":")]
            paths += [os.path.join(os.path.dirname(os.path.realpath(
                            inspect.getfile(inspect.currentframe()))),"limits")]

        # Load the limits.
        for path in (paths,) if isinstance(paths, str) else paths:
            limits = sorted(os.listdir(path))
            for limit in limits:
                if not limit.endswith(".py"): continue
                try: self[limit[0:-3]] = Limit(limit[0:-3], path)
                except: pass
