# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit corresponds to figure 3a (filled purple curve labeled CHARM)
of Tsai:2019mtm and was provided directly by the authors.

Production is from both eta and eta', but eta production is assumed to
provide the dominant contribution. This assumption can be tested by
modifying the relevant production fractions, which may also be made
mass dependent. Even if the eta' contributes to half the production
fraction, the change in limits is negligible.

This is a displaced search where the decay volume length over the
shielding length is 10/480.
"""
bibtex = """
@article{Tsai:2019mtm,
 author         = "Tsai, Yu-Dai and deNiverville, Patrick and Liu, Ming
                   Xiong",
 title          = "{The High-Energy Frontier of the Intensity Frontier:
                   Closing the Dark Photon, Inelastic Dark Matter, and Muon
                   g-2 Windows}",
 year           = "2019",
 eprint         = "1908.07525",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "FERMILAB-PUB-19-393-A-PPD",
 SLACcitation   = "%%CITATION = ARXIV:1908.07525;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production({"eta_gamma": 1.0, "eta'_gamma": 0.0})
decay = "e_e"
bounds = darkcast.Datasets("limits/CHARM_Tsai2019mtm.lmt")
efficiency = darkcast.Efficiency(lratio = 10.0/480.0)
