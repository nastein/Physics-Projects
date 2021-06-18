# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2018 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for LHCb searches using a di-electron final
state from D*0 -> D0 A' decays. This limit was extracted from figure 2
(blue lines) of Ilten:2015hya.

The limit is prompt, so t0 = 0. The selection is similar to the
phenomenology study of Ilten:2016tkc and so the definition of the
maximum proper lifetime from this paper is used, 0.04/(m - 2m_mu) +
0.1 ps.
"""
bibtex = """
@article{Ilten:2015hya,
 author         = "Ilten, Philip and Thaler, Jesse and Williams, Mike and
                   Xue, Wei",
 title          = "{Dark photons from charm mesons at LHCb}",
 journal        = "Phys. Rev.",
 volume         = "D92",
 year           = "2015",
 number         = "11",
 pages          = "115017",
 doi            = "10.1103/PhysRevD.92.115017",
 eprint         = "1509.06765",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "MIT-CTP-4702",
 SLACcitation   = "%%CITATION = ARXIV:1509.06765;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("D*0_D0")
decay = "e_e"
bounds = darkcast.Dataset("reach/LHCb_Ilten2015hya_prompt.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: (0.004/(m - 2*darkcast.pars.mfs['mu']) + 0.1)*1e-12)
