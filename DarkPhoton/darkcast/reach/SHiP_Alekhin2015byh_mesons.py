# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """ 
This limit is a projection for SHiP searches using meson decay
production. This limit was extracted from figure 2.6 (blue line) of
Alekhin:2015byh. Here the decay of the pi0 is taken as the predominant
production mechanism below the pi0 mass, and the decay of the eta
above.

This is a displaced search where the decay volume length over the
shielding length is 48/64.
"""
bibtex = """
@article{Alekhin:2015byh,
 author         = "Alekhin, Sergey and others",
 title          = "{A facility to Search for Hidden Particles at the CERN
                   SPS: the SHiP physics case}",
 journal        = "Rept. Prog. Phys.",
 volume         = "79",
 year           = "2016",
 number         = "12",
 pages          = "124201",
 doi            = "10.1088/0034-4885/79/12/124201",
 eprint         = "1504.04855",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "CERN-SPSC-2015-017, SPSC-P-350-ADD-1",
 SLACcitation   = "%%CITATION = ARXIV:1504.04855;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production({
    "pi0_gamma": lambda m: 1.0 if m <= darkcast.pars.mms["pi0"] else 0.0,
    "eta_gamma": lambda m: 1.0 if m > darkcast.pars.mms["pi0"] else 0.0})
production.name = "pi0_gamma"
decay = "e_e"
bounds = darkcast.Datasets("reach/SHiP_Alekhin2015byh_mesons.lmt")
efficiency = darkcast.Efficiency(lratio = 48.0/64.0)
