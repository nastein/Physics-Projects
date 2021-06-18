# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for the DarkLight search. A single prompt
bound is provided, 'DarkLight_Kahn2012br.lmt'. These limits are
extracted from figure 18 (right plot, purple line) of Kahn:2012br.

This search is prompt, and is not sensitive to X bosons with lifetimes
large enough to qualify as non-prompt; the efficiency ratio is assumed
to be unity, e.g. t0 = 0 and t1 = infinity.
"""
bibtex = """
@article{Kahn:2012br,
 author         = "Kahn, Yonatan and Thaler, Jesse",
 title          = "{Searching for an invisible A' vector boson with
                   DarkLight}",
 journal        = "Phys. Rev.",
 volume         = "D86",
 year           = "2012",
 pages          = "115012",
 doi            = "10.1103/PhysRevD.86.115012",
 eprint         = "1209.0777",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "MIT-CTP-4398",
 SLACcitation   = "%%CITATION = ARXIV:1209.0777;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Dataset("reach/DarkLight_Kahn2012br.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
