# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 2 (purple dotted curve labeled
E141) of Andreas:2012mt.

This is a displaced search where the decay volume length over the
shielding length is 35/0.12.
"""
bibtex = """
@article{Riordan:1987aw,
 author         = "Riordan, E. M. and others",
 title          = "{A Search for Short Lived Axions in an Electron Beam Dump
                   Experiment}",
 journal        = "Phys. Rev. Lett.",
 volume         = "59",
 year           = "1987",
 pages          = "755",
 doi            = "10.1103/PhysRevLett.59.755",
 reportNumber   = "SLAC-PUB-4280, UR-993, FERMILAB-PUB-87-251",
 SLACcitation   = "%%CITATION = PRLTA,59,755;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/E141_Riordan1987aw.lmt")
efficiency = darkcast.Efficiency(lratio = 35.0/0.12)
