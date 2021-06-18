# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import math
from . import utils, pars, model

###############################################################################
class BreitWignerError(Exception):
    """
    Simple exception for the 'BreitWigner' class.
    """
    pass

###############################################################################
class BreitWigner:
    """
    Provide a running width Breit-Wigner for a given resonance.
    
    mr:  mass of the resonance (GeV).
    wr:  width of the resonance (1/GeV).
    lr:  orbital angular momentum for the resonance.
    drs: decay channels for the resonance of the form 
         [branching, mass 0, mass 1].
    """
    ###########################################################################
    def __init__(self, resonance, lr = 1):
        """
        Initialize the needed data for a given resonance and orbital
        angular momentum.

        resonance: resonance string, e.g. 'rho0'.
        lr:        orbital angular momentum for the resonance, e.g. 0 for 
                   s-wave, 1 for p-wave, etc. If 'None', a fixed width is used.
        """
        self.lr = None if pars.bw == "fix" else lr
        try:
            self.mr, self.wr = pars.mms[resonance], pars.wms[resonance]
            self.sr, self.drs = self.mr**2, pars.dms[resonance]
        except:
            raise BreitWignerError(
                "No data is available for the %s resonance." % resonance)

    ###########################################################################
    def __call__(self, m):
        """
        Return the Breit-Wigner for a given mass.
        
        m: mass (GeV).
        """
        # Running width, equation A.3.
        if self.lr != None:
            k, sm = 0, m**2
            for br, m0, m1 in self.drs:
                if m > m0 + m1: k += br*math.sqrt(
                    ((sm - (m0 + m1)**2)*(sm - (m0 - m1)**2)/(4*sm))/
                    ((self.sr - (m0 + m1)**2)*(self.sr - (m0 - m1)**2)/
                     (4*self.sr)))**(2*self.lr + 1)
            return self.sr/(self.sr - sm - complex(0, self.sr*self.wr*k/m))
       
       # Fixed width, equation A.2.
        else:
            return self.sr/(self.sr - m**2 - complex(0, m*self.wr))

###############################################################################
class ProductionError(Exception):
    """
    Simple exception for the 'Production' class.
    """
    pass

###############################################################################
class Production:
    """
    Represents the possible production mechanisms for an X-boson,
    given an experiment.

    name:     name of the production.
    channels: list of production channels taking the form 
              [production of type 'Production', production fraction function]
    """
    ###########################################################################
    def __init__(self, channels, frac = 1.0):
        """
        Load a production, given its channel or channels. When
        specifying a single channel, 'channels' can be a mechanism
        name, e.g. 'p_brem' for proton-beam bremsstrahlung. The
        built-in mechanisms are:

        p_brem: proton-beam bremsstrahlung.
        A_brem: A-beam bremsstrahlung, where A can be any fundamental fermion.
        A_A:    Drell-Yan, where A can be any fundamental fermion.
        A:      vector meson mixing, where A is the vector meson.
        A_B:    meson decays of the form A -> B + X, where A is the decaying 
                meson, B is the SM daughter and X is the NP daughter.

        Alternatively, a user defined function that is mass and model
        dependent can be provided, taking the form 'mechanism(mass
        (GeV), model)'. The mechanism is assumed to be dependent upon
        the square of the global coupling.

        If multiple channels are specified, these are provided in a
        dictionary via 'channels' where the keys are the mechanisms
        and their associated values are the fractions. The mechanism
        is given as above. The fraction can be given either as a
        number, or a mass dependent function, e.g. 'fraction(mass
        (GeV))'.

        The 'Datasets' class provides a dictionary of datasets, and
        consequently can be passed as the 'channels' argument. The
        first column of the dataset is the mass interpolation points
        and all remaining columns are the mass dependent ratios for
        each mechanism. The first row specifies the built-in mechanism
        for each column. User defined mechanisms cannot be used here.
        """
        # Initialize the cached results.
        self.name = 'undefined'
        self.__cache = (None, None)

        # Multiple mechanisms from a dictionary.
        try:
            self.channels = []
            for prd, frc in channels.items():
                self.channels.append(Production(prd, frc))
            return
        except: pass

        # Pre-defined mechanism.
        if isinstance(channels, str):
            self.name = channels
            moms = channels.split('_')

            # Proton-beam bremsstrahlung, equation 2.4.
            if channels == "p_brem":
                self.__sigma = lambda m, model: (
                    2*model.xfs["u"](m) + model.xfs["d"](m))**2
    
            # Drell-Yan or lepton-beam bremsstrahlung, equation 2.3 and 2.6.
            elif (len(moms) == 2 and (moms[1] == 'brem' or moms[0] == moms[1]) 
                  and moms[0] in pars.mfs):
                self.__sigma = lambda m, model: (
                    model.xfs[moms[0]](m))**2
    
            # Vector meson decay, equation 2.11.
            elif channels in pars.rvs:
                self.__sigma = lambda m, model: (
                    model.trq(m, pars.tms[channels]))**2
    
            # Meson decay of the form A -> B + X, equation 2.7 - 2.10.
            elif len(moms) == 2 and moms[0] in pars.tms and moms[1] in pars.tms:
                ta, tb, vs = pars.tms[moms[0]], pars.tms[moms[1]], []
                for v in pars.rvs:
                    tv = pars.tms[v]
                    pf = utils.trace(ta, tb, tv)
                    if pf: vs.append((tv, pf, BreitWigner(v, 1)))
                self.__sigma = lambda m, model: (
                    abs(sum([pf*model.trq(m, tv)*
                             bw(m) for tv, pf, bw in vs])))**2
                    
            # Unknown mechanism.
            else: raise ProductionError(
                "Unknown production mechanism '%s'." % self.name)

        # User supplied mechanism.
        else:
            self.__sigma = channels
            self.__sigma(0, model.Model("dark_photon"))
                         
        # Set the channels.
        try: float(frac); self.__frac = lambda m, frac = frac: float(frac)
        except: self.__frac = frac
        self.channels = [self]

    ###########################################################################
    def ratio(self, m, g0, g1, model0, model1):
        """
        Return the cross-section ratio between two models, for a given
        mass and global couplings.
        
        m:      mass (GeV).
        g0:     global coupling for the first model.
        g1:     global coupling for the second model.
        model0: first model, numerator.
        model1: second model, denominator.
        """
        # Return the cached result if valid.
        if self.__cache[0] == m: return (g0/g1)**2*self.__cache[-1]

        # Calculate the result, equation 2.12.
        ratio = 0
        for channel in self.channels:
            den = channel.__sigma(m, model1)
            if den: ratio += channel.__frac(m)*channel.__sigma(m, model0)/den
        self.__cache = (m, ratio)
        return (g0/g1)**2*ratio
