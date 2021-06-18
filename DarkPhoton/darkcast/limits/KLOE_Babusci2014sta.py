# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 9 (blue filled curve labeled
KLOE(2)) of Anastasi:2016ktq.

This is a prompt search with a flight distance less than 8 cm. The
production is e+ e- -> gamma X at a center-of-mass energy of the phi
meson, so the gamma factor is on the order of (m_phi^2 + m^2)/(2 m_phi m).
"""
bibtex = """
@article{Babusci:2014sta,
 author         = "Babusci, D. and others",
 title          = "{Search for light vector boson production in $e^+e^-
                   \rightarrow \mu^+ \mu^- \gamma$ interactions with the KLOE
                   experiment}",
 collaboration  = "KLOE-2",
 journal        = "Phys. Lett.",
 volume         = "B736",
 year           = "2014",
 pages          = "459-464",
 doi            = "10.1016/j.physletb.2014.08.005",
 eprint         = "1404.7772",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1404.7772;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = "mu_mu"
bounds = darkcast.Dataset('limits/KLOE_Babusci2014sta.lmt')
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: 0.08/(
        (darkcast.pars.mms["phi"]**2 + m**2)
        /(2*darkcast.pars.mms["phi"]*m)*darkcast.pars.c))
