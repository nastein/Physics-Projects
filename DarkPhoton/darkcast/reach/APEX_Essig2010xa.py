# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for the APEX search. A single prompt bound
is provided, 'APEX_Essig2010xa.lmt'. These limits are extracted from
figure 1 (blue line) of Essig:2010xa.

No detailed information is given on the prompt-like
requirements. However, since the same coupling is used in production
and decay, the efficiency ratio is assumed to be unity, e.g. t0 = 0
and t1 = infinity.
"""
bibtex = """
@article{Essig:2010xa,
 author         = "Essig, Rouven and Schuster, Philip and Toro, Natalia and
                   Wojtsekhowski, Bogdan",
 title          = "{An Electron Fixed Target Experiment to Search for a New
                   Vector Boson A' Decaying to e+e-}",
 journal        = "JHEP",
 volume         = "02",
 year           = "2011",
 pages          = "009",
 doi            = "10.1007/JHEP02(2011)009",
 eprint         = "1001.2557",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "SLAC-PUB-13882, SU-ITP-10-01",
 SLACcitation   = "%%CITATION = ARXIV:1001.2557;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Dataset("reach/APEX_Essig2010xa.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
