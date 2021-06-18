# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (solid gray curve labeled NOMAD
& PS191 just below the nu-Cal I (pi0) line) of Blumlein:2013cua.

This is a displaced search where the decay volume length over the
shielding length is 7.5/835.
"""
bibtex = """
@article{Astier:2001ck,
 author         = "Astier, P. and others",
 title          = "{Search for heavy neutrinos mixing with tau neutrinos}",
 collaboration  = "NOMAD",
 journal        = "Phys. Lett.",
 volume         = "B506",
 year           = "2001",
 pages          = "27-38",
 doi            = "10.1016/S0370-2693(01)00362-8",
 eprint         = "hep-ex/0101041",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-EP-2001-005",
 SLACcitation   = "%%CITATION = HEP-EX/0101041;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("pi0_gamma")
decay = "e_e"
bounds = darkcast.Datasets("limits/NOMAD_Astier2001ck.lmt")
efficiency = darkcast.Efficiency(lratio = 7.5/835.0)
