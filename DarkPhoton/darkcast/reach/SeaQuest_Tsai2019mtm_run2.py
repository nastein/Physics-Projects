# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This projection corresponds to figure 3b (dashed green curve labeled
SQ/DQ (1e20 POT)) of Tsai:2019mtm and was provided directly by the
authors. A total of 1e20 protons on target is assumed for run 2 data
taking.

Production is from both pi0 and eta meson decays. The production
fraction is assumed to be entirely from pi0 decays when kinematically
available, otherwise eta decays.

In Tsai:2019mtm this is considered a displaced search with a decay
volume of 7 meters and shielding length of 5 meters.
"""
bibtex = """
@article{Gardner:2015wea,
 author         = "Gardner, S. and Holt, R. J. and Tadepalli, A. S.",
 title          = "{New Prospects in Fixed Target Searches for Dark Forces
                   with the SeaQuest Experiment at Fermilab}",
 journal        = "Phys. Rev.",
 volume         = "D93",
 year           = "2016",
 number         = "11",
 pages          = "115015",
 doi            = "10.1103/PhysRevD.93.115015",
 eprint         = "1509.00050",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 SLACcitation   = "%%CITATION = ARXIV:1509.00050;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production({
    "pi0_gamma": lambda m: 1.0*(m <  0.13498 - 2*5.110e-04),
    "eta_gamma": lambda m: 1.0*(m >= 0.13498 - 2*5.110e-04)})
decay = ["e_e", "mu_mu"]
bounds = darkcast.Datasets("reach/SeaQuest_Tsai2019mtm_run2.lmt")
efficiency = darkcast.Efficiency(lratio = 7.0/5.0)
