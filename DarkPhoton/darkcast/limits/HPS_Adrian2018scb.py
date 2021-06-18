# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 3 (blue curve labeled
2015 engineering run) of Adrian:2018scb.

No detailed information is given on the prompt-like
requirements. However, since the same coupling is used in production
and decay, the efficiency ratio is assumed to be unity, e.g. t0 = 0
and t1 = infinity.
"""
bibtex = """
@article{Adrian:2018scb,
 author         = "Adrian, P. H. and others",
 title          = "{Search for a Dark Photon in Electro-Produced
                   $e^{+}e^{-}$ Pairs with the Heavy Photon Search Experiment
                   at JLab}",
 year           = "2018",
 eprint         = "1807.11530",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1807.11530;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Dataset("limits/HPS_Adrian2018scb.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
