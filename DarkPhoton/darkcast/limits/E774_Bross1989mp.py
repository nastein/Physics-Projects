# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 2 (orange dashed curve labeled
E774) of Andreas:2012mt.

This is a displaced search where the decay volume length over the
shielding length is 2/0.3.
"""
bibtex = """
@article{Bross:1989mp,
 author         = "Bross, A. and Crisler, M. and Pordes, Stephen H. and
                   Volk, J. and Errede, S. and Wrbanek, J.",
 title          = "{A Search for Short-lived Particles Produced in an
                   Electron Beam Dump}",
 journal        = "Phys. Rev. Lett.",
 volume         = "67",
 year           = "1991",
 pages          = "2942-2945",
 doi            = "10.1103/PhysRevLett.67.2942",
 reportNumber   = "FERMILAB-PUB-89-138-E",
 SLACcitation   = "%%CITATION = PRLTA,67,2942;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/E774_Bross1989mp.lmt")
efficiency = darkcast.Efficiency(lratio = 2.0/0.3)
