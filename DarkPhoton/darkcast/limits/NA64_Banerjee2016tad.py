# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (blue curve labeled NA64) of
Lees:2017lec.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{Banerjee:2016tad,
 author         = "Banerjee, D. and others",
 title          = "{Search for invisible decays of sub-GeV dark photons in
                   missing-energy events at the CERN SPS}",
 collaboration  = "NA64",
 journal        = "Phys. Rev. Lett.",
 volume         = "118",
 year           = "2017",
 number         = "1",
 pages          = "011802",
 doi            = "10.1103/PhysRevLett.118.011802",
 eprint         = "1610.02988",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1610.02988;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                           99.0*model.width("visible", m))
production = darkcast.Production("e_brem")
decay = "invisible"
bounds = darkcast.Dataset("limits/NA64_Banerjee2016tad.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
