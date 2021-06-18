# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 2 (green dash-dotted curve labeled
KEK) of Andreas:2012mt.

This is a displaced search where the decay volume length over the
shielding length is 2.2/2.4.
"""
bibtex = """
@article{Konaka:1986cb,
 author         = "Konaka, A. and others",
 title          = "{Search for Neutral Particles in Electron Beam Dump
                   Experiment}",
 booktitle      = "{Proceedings, 23RD International Conference on High
                   Energy Physics, JULY 16-23, 1986, Berkeley, CA}",
 journal        = "Phys. Rev. Lett.",
 volume         = "57",
 year           = "1986",
 pages          = "659",
 doi            = "10.1103/PhysRevLett.57.659",
 reportNumber   = "KEK-Preprint-86-9",
 SLACcitation   = "%%CITATION = PRLTA,57,659;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/KEK_Konaka1986cb.lmt")
efficiency = darkcast.Efficiency(lratio = 2.2/2.4)
