# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (light blue curve labeled
NA64) of Banerjee:2019hmi.

This is a displaced search where the decay volume length over the
shielding length is approximately 4/1.
"""
bibtex = """
@article{Banerjee:2019hmi,
 author         = "Banerjee, D. and others",
 title          = "{Improved limits on a hypothetical X(16.7) boson and a
                   dark photon decaying into $e^+e^-$ pairs}",
 year           = "2019",
 eprint         = "1912.11389",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-EP-2019-284",
 SLACcitation   = "%%CITATION = ARXIV:1912.11389;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/NA64_Banerjee2019hmi.lmt")
efficiency = darkcast.Efficiency(lratio = 4.0/1.0)
