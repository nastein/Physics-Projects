# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for Belle II searches. A single prompt
bound is provided, 'BelleII_Kou2018nap.lmt'. These limits were
provided by the authors although additional limits can be found in
figure 201 of Kou:2018nap.

This search is prompt, and is not sensitive to X bosons with lifetimes
large enough to qualify as non-prompt; the efficiency ratio is assumed
to be unity, e.g. t0 = 0 and t1 = infinity.
"""
bibtex = """
@article{Kou:2018nap,
 author         = "Kou, E. and others",
 title          = "{The Belle II Physics Book}",
 year           = "2018",
 eprint         = "1808.10567",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "KEK Preprint 2018-27, BELLE2-PUB-PH-2018-001,
                   FERMILAB-PUB-18-398-T, JLAB-THY-18-2780, INT-PUB-18-047,
                   UWThPh 2018-26",
 SLACcitation   = "%%CITATION = ARXIV:1808.10567;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = ["e_e", "mu_mu"]
bounds = darkcast.Dataset("reach/BelleII_Kou2018nap.lmt")
efficiency = darkcast.Efficiency(t0 = 0, t1 = float("inf"))
