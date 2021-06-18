# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 3 (light blue curve) of
Adachi:2019otg, divided by a factor of 'pars.ge' given by 3.02822e-1
to match the Darkcast dark photon model.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{Adachi:2019otg,
 author         = "Adachi, I. and others",
 title          = "{Search for an invisibly decaying $Z^{\prime}$ boson at
                   Belle II in $e^+ e^- \to \mu^+ \mu^- (e^{\pm} \mu^{\mp})$
                   plus missing energy final states}",
 collaboration  = "Belle II",
 year           = "2019",
 eprint         = "1912.11276",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1912.11276;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                       99.0*model.width("visible", m))
production = darkcast.Production("mu_brem")
decay = "invisible"
bounds = darkcast.Dataset("limits/BelleII_Adachi2019otg.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
