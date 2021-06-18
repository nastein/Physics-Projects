# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for the MESA search. A single prompt bound
is provided, 'MESA_Beranek2013yqa.lmt'. These limits are extracted
from figure 12 (blue and red lines) of Beranek:2013yqa.

This search is prompt, and is not sensitive to X bosons with lifetimes
large enough to qualify as non-prompt; the efficiency ratio is assumed
to be unity, e.g. t0 = 0 and t1 = infinity.
"""
bibtex = """
@article{Beranek:2013yqa,
 author         = "Beranek, T. and Merkel, H. and Vanderhaeghen, M.",
 title          = "{Theoretical framework to analyze searches for hidden
                   light gauge bosons in electron scattering fixed target
                   experiments}",
 journal        = "Phys. Rev.",
 volume         = "D88",
 year           = "2013",
 pages          = "015032",
 doi            = "10.1103/PhysRevD.88.015032",
 eprint         = "1303.2540",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 SLACcitation   = "%%CITATION = ARXIV:1303.2540;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = ["e_e", "mu_mu"]
bounds = darkcast.Dataset("reach/MESA_Beranek2013yqa.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
