# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 3 (blue curve labeled NA64) of
NA64:2019imj.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{NA64:2019imj,
 author         = "Banerjee, D. and others",
 title          = "{Dark matter search in missing energy events with NA64}",
 journal        = "Phys. Rev. Lett.",
 volume         = "123",
 year           = "2019",
 number         = "12",
 pages          = "121801",
 doi            = "10.1103/PhysRevLett.123.121801",
 eprint         = "1906.00176",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-EP-2019-116",
 SLACcitation   = "%%CITATION = ARXIV:1906.00176;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                       99.0*model.width("visible", m))
production = darkcast.Production("e_brem")
decay = "invisible"
bounds = darkcast.Dataset("limits/NA64_NA642019imj.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
