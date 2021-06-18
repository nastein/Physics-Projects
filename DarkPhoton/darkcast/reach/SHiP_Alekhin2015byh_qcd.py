# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """ 
This limit is a projection for SHiP searches using QCD production, q
qbar -> V and q g -> q V. This limit was extracted from figure 2.6
(orange line) of Alekhin:2015byh. The flavor dependent production is
taken from the Drell-Yan component of LHC production at low mass from
Ilten:2016tk. This is a very rough approximation.

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
production = darkcast.Production({"u_u": 0.18, "d_d": 0.81, "s_s": 0.01})
production.name = "QCD"
decay = "e_e"
bounds = darkcast.Datasets("reach/SHiP_Alekhin2015byh_qcd.lmt")
efficiency = darkcast.Efficiency(lratio = 48.0/64.0)
