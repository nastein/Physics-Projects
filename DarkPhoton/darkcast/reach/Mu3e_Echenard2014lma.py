# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for the Mu3e search. A single prompt bound
is provided, 'Mu3e_Echenard2014lma.lmt'. These limits are extracted
from figure 7 (blue and red lines) of Echenard:2014lma.

The production is from muon decays and is approximated as production
from electron bremsstrahlung. However, there is also a contribution
from muon bremsstrahlung and W boson bremsstrahlung.

This search is prompt, and is not sensitive to X bosons with lifetimes
large enough to qualify as non-prompt; the efficiency ratio is assumed
to be unity, e.g. t0 = 0 and t1 = infinity.
"""
bibtex = """
@article{Echenard:2014lma,
 author         = "Echenard, Bertrand and Essig, Rouven and Zhong, Yi-Ming",
 title          = "{Projections for Dark Photon Searches at Mu3e}",
 journal        = "JHEP",
 volume         = "01",
 year           = "2015",
 pages          = "113",
 doi            = "10.1007/JHEP01(2015)113",
 eprint         = "1411.1770",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "YITP-SB-14-36",
 SLACcitation   = "%%CITATION = ARXIV:1411.1770;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("mu_brem")
decay = "e_e"
bounds = darkcast.Dataset("reach/Mu3e_Echenard2014lma.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
