# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """ 
This limit is a projection for FASER searches using meson decay
production. This limit was extracted from the right plot of figure 6
(black line labeled FASER 2) of Ariga:2018uku. Because this limit
extends well above the mass of the eta meson it is necessary to
account for production from sources other than pi0 and eta
decays. Consequently, the production ratios from Ilten:2016tkc are
used. While these were specifically calculated for LHCb acceptance
they should also provide a decent approximation for FASER acceptance.

This is a displaced search where the decay volume length over the
shielding length is 5/480.
"""
bibtex = """
@article{Ariga:2018uku,
 author         = "Ariga, Akitaka and others",
 title          = "{FASER's Physics Reach for Long-Lived Particles}",
 collaboration  = "FASER",
 year           = "2018",
 eprint         = "1811.12522",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "UCI-TR-2018-19, KYUSHU-RCAPP-2018-06",
 SLACcitation   = "%%CITATION = ARXIV:1811.12522;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production(
    darkcast.Datasets("limits/LHCb_Aaij2017rft.prd"))
decay = "e_e"
bounds = darkcast.Datasets("reach/FASER_Ariga2018uku_hllhc.lmt")
efficiency = darkcast.Efficiency(lratio = 5.0/480.0)
