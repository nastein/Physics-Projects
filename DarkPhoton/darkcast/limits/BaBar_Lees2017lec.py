# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (green curve labeled BaBar) of
Lees:2017lec.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{Lees:2017lec,
 author         = "Lees, J. P. and others",
 title          = "{Search for Invisible Decays of a Dark Photon Produced in
                   ${e}^{+}{e}^{-}$ Collisions at BaBar}",
 collaboration  = "BaBar",
 journal        = "Phys. Rev. Lett.",
 volume         = "119",
 year           = "2017",
 number         = "13",
 pages          = "131804",
 doi            = "10.1103/PhysRevLett.119.131804",
 eprint         = "1702.03327",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "BABAR-PUB-17-001, SLAC-PUB-16923",
 SLACcitation   = "%%CITATION = ARXIV:1702.03327;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                       99.0*model.width("visible", m))
production = darkcast.Production("e_e")
decay = "invisible"
bounds = darkcast.Dataset("limits/BaBar_Lees2017lec.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
