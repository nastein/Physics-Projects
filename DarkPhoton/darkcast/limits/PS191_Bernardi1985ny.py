# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (tip of the solid gray curve
labeled NOMAD & PS191) of Blumlein:2013cua.

This is a displaced search where the decay volume length over the
shielding length is 7/128.
"""
bibtex = """
@article{Bernardi:1985ny,
      author         = "Bernardi, G. and others",
      title          = "{Search for Neutrino Decay}",
      journal        = "Phys. Lett.",
      volume         = "166B",
      year           = "1986",
      pages          = "479-483",
      doi            = "10.1016/0370-2693(86)91602-3",
      reportNumber   = "CERN-EP/85-177",
      SLACcitation   = "%%CITATION = PHLTA,166B,479;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("pi0_gamma")
decay = "e_e"
bounds = darkcast.Datasets("limits/PS191_Bernardi1985ny.lmt")
efficiency = darkcast.Efficiency(lratio = 7.0/128.0)
