# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 9 (red filled curve labeled
KLOE(4)) of Anastasi:2016ktq.

This is a prompt search with a flight distance less than 8 cm. The
production is e+ e- -> gamma X at a center-of-mass energy of the phi
meson, so the gamma factor is on the order of (m_phi^2 + m^2)/(2 m_phi m).
"""
bibtex = """
@article{Anastasi:2016ktq,
 author         = "Anastasi, A. and others",
 title          = "{Limit on the production of a new vector boson in
                   $\mathrm{e^+ e^-}\rightarrow {\rm U}\gamma$, U$\rightarrow
                   \pi^+\pi^-$ with the KLOE experiment}",
 collaboration  = "KLOE-2",
 journal        = "Phys. Lett.",
 volume         = "B757",
 year           = "2016",
 pages          = "356-361",
 doi            = "10.1016/j.physletb.2016.04.019",
 eprint         = "1603.06086",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1603.06086;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = "pi+_pi-"
bounds = darkcast.Dataset("limits/KLOE_Anastasi2016ktq.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: 0.08/(
        (darkcast.pars.mms["phi"]**2 + m**2)
        /(2*darkcast.pars.mms["phi"]*m)*darkcast.pars.c))
