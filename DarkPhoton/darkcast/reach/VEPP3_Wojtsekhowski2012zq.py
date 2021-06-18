# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for the VEPP-3 search for dark photon
produced from a positron beam incident on a hydrogen target. This
limit is extracted from figure 9 (magenta line) of
Wojtsekhowski:2012zq.

No detailed information is given on the prompt-like
requirements. However, since the same coupling is used in production
and decay, the efficiency ratio is assumed to be unity, e.g. t0 = 0
and t1 = infinity.
"""
bibtex = """
@article{Wojtsekhowski:2012zq,
 author         = "Wojtsekhowski, B. and Nikolenko, D. and Rachek, I.",
 title          = "{Searching for a new force at VEPP-3}",
 year           = "2012",
 eprint         = "1207.5089",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "JLAB-PHY-12-1597, DOE-OR-23177-2213",
 SLACcitation   = "%%CITATION = ARXIV:1207.5089;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = "e_e"
bounds = darkcast.Dataset("reach/VEPP3_Wojtsekhowski2012zq.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
