# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 3 (blue curve labeled
NA64) of Banerjee:2018vgk.

This is a displaced search where the decay volume length over the
shielding length is approximately 4/1.
"""
bibtex = """
@article{Banerjee:2018vgk,
 author         = "Banerjee, D. and others",
 title          = "{Search for a new X(16.7) boson and dark photons in the
                   NA64 experiment at CERN}",
 collaboration  = "NA64",
 year           = "2018",
 eprint         = "1803.07748",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-EP-2018-043, CERN-EP-2018-043",
 SLACcitation   = "%%CITATION = ARXIV:1803.07748;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/NA64_Banerjee2018vgk.lmt")
efficiency = darkcast.Efficiency(lratio = 4.0/1.0)
