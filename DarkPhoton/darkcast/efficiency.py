# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import math
from . import utils

###############################################################################
class EfficiencyError(Exception):
    """
    Simple exception for the 'Efficiency' class.
    """
    pass

###############################################################################
class Efficiency:
    """
    Provides the efficiency for a given experiment.
    """
    ###########################################################################
    def __init__(self, t0 = 0, t1 = None, lratio = None, rvals = False):
        """
        Initialize an efficiency. The efficiency is defined in terms
        of a proper lifetime interval between t0 and t1 (seconds). For
        a prompt experiment, only the upper proper time needs to be
        provided, e.g. t1. For a displaced experiment, only 'lratio'
        (e.g. the ratio between the decay volume length and the shielding
        length, L_det/L_sh) needs to be provided, and t0 and t1 will
        be solved for using the lower and upper limits.

        When r-value limits are given, i.e. the ratio between observed
        and expected as a function of mass and lifetime, 'rvals'
        should be set as true. Here, the efficiency ratio for
        lifetimes above the maximum is approximated by t0 = 0 and t1 =
        tau_max, while below the minimum by t0 = tau_min and t1 =
        infinity.

        t0:     value or function for the lower proper lifetime (seconds).
        t1:     value or function for the upper proper lifetime (seconds).
        lratio: ratio between the decay and shielding volume for a beam-dump
                experiment.
        rvals:  true if the efficiency is for limits with r-values, false
                otherwise.
        """
        # Return if efficiency for a r-value limit.
        self.__rvals = rvals
        if self.__rvals: return

        # Initialize the cached results.
        self.__cache = (None, None, None)

        # Set the lower proper time.
        if lratio != None: self.__lratio = 1.0 + lratio
        else:
            self.__lratio = None
            try: t0 = float(t0); self.__t0 = lambda m: t0
            except: self.__t0 = t0
 
        # Set the upper proper time.
        if lratio == None:
            try: t1 = float(t1); self.__t1 = lambda m: t1
            except: self.__t1 = t1

    ###########################################################################
    def __ts(self, m, limit):
        """
        Return the proper time fiducial. If a displaced limit and
        'lratio' was specified, solve for t0 and t1.

        m:     mass (GeV).
        limit: 'Limit' which includes a model and lower/upper bounds.
        """
        # Cached proper times.
        if self.__cache[0:-1] == (m, limit): return self.__cache[-1]

        # Fiducial from displaced limits with no shielding.
        if self.__lratio == float("inf"):
            t0 = limit.model.tau(m, limit.bounds["lower"](m))
            t1 = limit.model.tau(m, limit.bounds["upper"](m))
        
        # Fiducial from displaced limit.
        elif self.__lratio:

            # Solve t0 from equation 2.24.
            g0 = limit.bounds["lower"](m)
            g1 = limit.bounds["upper"](m)
            tau0 = limit.model.tau(m, g0)
            tau1 = limit.model.tau(m, g1)
            f = lambda t: (
                g1**2*(math.exp(-t/tau1) - math.exp(-t*self.__lratio/tau1))
                - g0**2*(math.exp(-t/tau0) - math.exp(-t*self.__lratio/tau0)))
            t0 = utils.solve(f)

            # Set t1 from equation 2.22.
            t1 = t0*self.__lratio

        # Fiducial from user defined proper times.
        else: t0, t1 = self.__t0(m), self.__t1(m)
        self.__cache = (m, limit, (t0, t1))
        return t0, t1
    
    ###########################################################################
    def ratio(self, m, limit, tau0, tau1):
        """
        Return the efficiency ratios for a given mass, limit, and
        lifetimes.
        
        m:     mass (GeV).
        limit: 'Limit' which includes a model and lower/upper bounds.
        tau0:  numerator lifetime (seconds).
        tau1:  denominator lifetime (seconds).
        """
        # If r-value limit, use equation C.4.
        if self.__rvals:
            t0, t1 = (0, tau1) if tau0 > tau1 else (tau1, float("inf"))
        else: t0, t1 = self.__ts(m, limit)

        # Ratio of efficiencies, given by equation 2.23.
        num = math.exp(-t0/tau0) - math.exp(-t1/tau0)
        den = math.exp(-t0/tau1) - math.exp(-t1/tau1)
        return num/den if den else float("inf")
