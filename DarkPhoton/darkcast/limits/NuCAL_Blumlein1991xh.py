# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 5 (solid red curve labeled
nu-Cal I (p-Bremstrahlung)) of Blumlein:2013cua.

This is a displaced search where the decay volume length over the
shielding length is 23/64.
"""
bibtex = """
@article{Blumlein:1991xh,
 author         = "Blumlein, J. and others",
 title          = "{Limits on the mass of light (pseudo)scalar particles
                   from Bethe-Heitler e+ e- and mu+ mu- pair production in a
                   proton - iron beam dump experiment}",
 journal        = "Int. J. Mod. Phys.",
 volume         = "A7",
 year           = "1992",
 pages          = "3835-3850",
 doi            = "10.1142/S0217751X9200171X",
 reportNumber   = "PHE-91-11",
 SLACcitation   = "%%CITATION = IMPAE,A7,3835;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("p_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/NuCAL_Blumlein1991xh.lmt")
efficiency = darkcast.Efficiency(lratio = 23.0/64.0)
