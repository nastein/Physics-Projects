# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 7 (red filled curve labeled
NA62) of the associated paper, CortinaGil:2019nuo.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{CortinaGil:2019nuo,
 author         = "Cortina Gil, Eduardo and others",
 title          = "{Search for production of an invisible dark photon in
                   $\pi^0$ decays}",
 collaboration  = "NA62",
 year           = "2019",
 eprint         = "1903.08767",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-EP-2019-048",
 SLACcitation   = "%%CITATION = ARXIV:1903.08767;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                           99.0*model.width("visible", m))
production = darkcast.Production("pi0_gamma")
decay = "invisible"
bounds = darkcast.Dataset("limits/NA62_CortinaGil2019nuo.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
