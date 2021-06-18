# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 8 (black dashed curve labeled
KLOE(mumu)) of Anastasi:2018azp.

This is a prompt search with a flight distance less than 8 cm. The
production is e+ e- -> gamma X at a center-of-mass energy of the phi
meson, so the gamma factor is on the order of (m_phi^2 + m^2)/(2 m_phi m).
"""
bibtex = """
@article{Anastasi:2018azp,
 author         = "Anastasi, A. and others",
 title          = "{Combined limit on the production of a light gauge boson
                   decaying into $\mu^+\mu^-$ and $\pi^+\pi^-$}",
 collaboration  = "KLOE-2",
 journal        = "Submitted to: Phys. Lett. B",
 year           = "2018",
 eprint         = "1807.02691",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 SLACcitation   = "%%CITATION = ARXIV:1807.02691;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_e")
decay = "mu_mu"
bounds = darkcast.Dataset('limits/KLOE_Anastasi2018azp_mumu.lmt')
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: 0.08/(
        (darkcast.pars.mms["phi"]**2 + m**2)
        /(2*darkcast.pars.mms["phi"]*m)*darkcast.pars.c))
