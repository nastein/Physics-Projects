# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (solid blue curve labeled
nu-Cal I (pi0)) of Blumlein:2013cua.

This is a displaced search where the decay volume length over the
shielding length is 23/64.
"""
bibtex = """
@article{Blumlein:1990ay,
 author         = "Blumlein, J. and others",
 title          = "{Limits on neutral light scalar and pseudoscalar
                   particles in a proton beam dump experiment}",
 journal        = "Z. Phys.",
 volume         = "C51",
 year           = "1991",
 pages          = "341-350",
 doi            = "10.1007/BF01548556",
 reportNumber   = "PHE-90-03",
 SLACcitation   = "%%CITATION = ZEPYA,C51,341;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("pi0_gamma")
decay = "e_e"
bounds = darkcast.Datasets("limits/NuCAL_Blumlein1990ay.lmt")
efficiency = darkcast.Efficiency(lratio = 23.0/64.0)
