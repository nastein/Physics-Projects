# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 15 (blue curve labeled NA64) of
Banerjee:2017hhz.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{Banerjee:2017hhz,
 author         = "Banerjee, D. and others",
 title          = "{Search for vector mediator of Dark Matter production in
                   invisible decay mode}",
 collaboration  = "NA64",
 journal        = "Phys. Rev.",
 volume         = "D97",
 year           = "2018",
 number         = "7",
 pages          = "072002",
 doi            = "10.1103/PhysRevD.97.072002",
 eprint         = "1710.00971",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1710.00971;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                       99.0*model.width("visible", m))
production = darkcast.Production("e_brem")
decay = "invisible"
bounds = darkcast.Dataset("limits/NA64_Banerjee2017hhz.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
