# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is from figure 8 (blue curve labeled LEP) of Ilten which is
recast from Fox:2011fx based on LEP data. The dark photon is assumed
to be on-shell, the dark matter coupling much larger than the
electromagnetic coupling, and the dark matter mass much less than the
dark photon mass.

This search is an invisible search and so the efficiency ratio is
assumed to be unity, e.g. t0 = 0 and t1 = infinity. Consequently, the
invisible width does not matter as long as it is proportional to the
visible width. Here, the invisible width is set as 99 times the
visible width, e.g. a branching fraction of 99%.
"""
bibtex = """
@article{Fox:2011fx,
 author         = "Fox, Patrick J. and Harnik, Roni and Kopp, Joachim and
                   Tsai, Yuhsin",
 title          = "{LEP Shines Light on Dark Matter}",
 journal        = "Phys. Rev.",
 volume         = "D84",
 year           = "2011",
 pages          = "014028",
 doi            = "10.1103/PhysRevD.84.014028",
 eprint         = "1103.0240",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "FERMILAB-PUB-11-039-T",
 SLACcitation   = "%%CITATION = ARXIV:1103.0240;%%"
}
"""
model = darkcast.Model("dark_photon", iwidth = lambda m, model:
                           99.0*model.width("visible", m))
production = darkcast.Production("e_e")
decay = "invisible"
bounds = darkcast.Dataset("limits/LEP_Fox2011fx.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
