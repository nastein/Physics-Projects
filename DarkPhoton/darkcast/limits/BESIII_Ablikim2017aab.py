# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 3 (yellow fill labeled BESIII) of
Ablikim:2017aab. In this limit both the electron and muon final states
were combined.

This search is prompt, and is not sensitive to X bosons with lifetimes
large enough to qualify as non-prompt; the efficiency ratio is assumed
to be unity, e.g. t0 = 0 and t1 = infinity.
"""
bibtex = """
@article{Ablikim:2017aab,
 author         = "Ablikim, M. and others",
 title          = "{Dark Photon Search in the Mass Range Between 1.5 and 3.4
                   GeV/$c^2$}",
 collaboration  = "BESIII",
 journal        = "Phys. Lett.",
 volume         = "B774",
 year           = "2017",
 pages          = "252-257",
 doi            = "10.1016/j.physletb.2017.09.067",
 eprint         = "1705.04265",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1705.04265;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = ["e_e", "mu_mu"]
bounds = darkcast.Dataset("limits/BESIII_Ablikim2017aab.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))

