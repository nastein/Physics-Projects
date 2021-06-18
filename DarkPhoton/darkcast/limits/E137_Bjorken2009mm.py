# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 1 (left, dark fill labeled
E137) of Bjorken:2009mm. The same limit, but determined using different
assumptions is available in E137_Andreas2012mt.

This is a displaced search where the decay volume length over the
shielding length is 204/179.
"""
bibtex = """
@article{Bjorken:2009mm,
 author         = "Bjorken, James D. and Essig, Rouven and Schuster, Philip 
                   and Toro, Natalia",
 title          = "{New Fixed-Target Experiments to Search for Dark Gauge 
                   Forces}",
 eprint         = "0906.0580",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "SLAC-PUB-13650, SU-ITP-09-22",
 doi            = "10.1103/PhysRevD.80.075018",
 journal        = "Phys. Rev. D",
 volume         = "80",
 pages          = "075018",
 year           = "2009"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/E137_Bjorken2009mm.lmt")
efficiency = darkcast.Efficiency(lratio = 204.0/179.0)
