# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 9 (blue filled curve labeled
KLOE(1)) of Anastasi:2016ktq. The dark photon is assumed to decay only
into e_e and mu_mu final states and so the model is modified accordingly.

This is a prompt search with a flight distance less than 8 cm. The
production is phi -> eta X, where the phi meson is produced at rest so 
the gamma factor is on the order of (m_phi^2 + m^2 - m_eta^2)/(2 m_phi m).
"""
bibtex = """
@article{Babusci:2012cr,
 author         = "Babusci, D. and others",
 title          = "{Limit on the production of a light vector gauge boson in
                   phi meson decays with the KLOE detector}",
 collaboration  = "KLOE-2",
 journal        = "Phys. Lett.",
 volume         = "B720",
 year           = "2013",
 pages          = "111-115",
 doi            = "10.1016/j.physletb.2013.01.067",
 eprint         = "1210.3927",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1210.3927;%%"
}
"""
model = darkcast.Model("dark_photon", states = ["e_e", "mu_mu"])
production = darkcast.Production("phi_eta")
decay = "e_e"
bounds = darkcast.Dataset("limits/KLOE_Babusci2012cr.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: 0.08/(
        (darkcast.pars.mms["phi"]**2 + m**2 - darkcast.pars.mms["eta"]**2)
        /(2*darkcast.pars.mms["phi"]*m)*darkcast.pars.c))
