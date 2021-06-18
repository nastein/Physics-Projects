# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (gray curve labeled
APEX) of Blumlein:2013cua.

No detailed information is given on the prompt-like
requirements. However, since the same coupling is used in production
and decay, the efficiency ratio is assumed to be unity, e.g. t0 = 0
and t1 = infinity.
"""
bibtex = """
@article{Abrahamyan:2011gv,
 author         = "Abrahamyan, S. and others",
 title          = "{Search for a New Gauge Boson in Electron-Nucleus
                   Fixed-Target Scattering by the APEX Experiment}",
 collaboration  = "APEX",
 journal        = "Phys. Rev. Lett.",
 volume         = "107",
 year           = "2011",
 pages          = "191804",
 doi            = "10.1103/PhysRevLett.107.191804",
 eprint         = "1108.2750",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "JLAB-PHY-11-1406, SLAC-PUB-14491,
                   JLAB-PHY-11-1406---SLAC-PUB-14491",
 SLACcitation   = "%%CITATION = ARXIV:1108.2750;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Dataset("limits/APEX_Abrahamyan2011gv.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
