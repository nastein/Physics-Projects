# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit corresponds to figure 3a (solid red curve labeled NA62)
of Tsai:2019mtm and was provided directly by the authors.

Production is from both pi0 and eta meson decays. The production
fraction is assumed to be entirely from pi0 decays when kinematically
available, otherwise eta decays.

NA62 has 82 meters of shielding followed by a decay volume of 135
meters. However, the first spectrometer chamber is 180 meters from the
interaction point, reducing the effective decay volume. Following the
convention of Tsai:2019mtm, a decay volume of 75 meters is used with a
shielding length of 142 meters.
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
production = darkcast.Production({
    "pi0_gamma": lambda m: 1.0*(m <  0.13498 - 2*5.110e-04),
    "eta_gamma": lambda m: 1.0*(m >= 0.13498 - 2*5.110e-04)})
decay = ["e_e", "mu_mu"]
bounds = darkcast.Datasets("reach/NA62_Tsai2019mtm.lmt")
efficiency = darkcast.Efficiency(lratio = 75.0/142.0)
