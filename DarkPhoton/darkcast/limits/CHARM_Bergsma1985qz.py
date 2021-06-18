# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (solid gray curve labeled
CHARM) of Blumlein:2013cua.

Production is from both eta and eta', but eta production is assumed to
provide the dominant contribution. This assumption can be tested by
modifying the relevant production fractions, which may also be made
mass dependent. Even if the eta' contributes to half the production
fraction, the change in limits is negligible.

This is a displaced search where the decay volume length over the
shielding length is 10/480.
"""
bibtex = """
@article{Bergsma:1985qz,
 author         = "Bergsma, F. and others",
 title          = "{Search for Axion Like Particle Production in 400-{GeV}
                   Proton - Copper Interactions}",
 collaboration  = "CHARM",
 journal        = "Phys. Lett.",
 volume         = "157B",
 year           = "1985",
 pages          = "458-462",
 doi            = "10.1016/0370-2693(85)90400-9",
 reportNumber   = "CERN-EP-85-38",
 SLACcitation   = "%%CITATION = PHLTA,157B,458;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production({"eta_gamma": 1.0, "eta'_gamma": 0.0})
decay = "e_e"
bounds = darkcast.Datasets("limits/CHARM_Bergsma1985qz.lmt")
efficiency = darkcast.Efficiency(lratio = 10.0/480.0)
