# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was extracted from figure 2 (red dashed curve labeled
E137) of Andreas:2012mt. The same limit, but determined using different
assumptions is available in E137_Bjorken:2009mm.

This is a displaced search where the decay volume length over the
shielding length is 204/179.
"""
bibtex = """
@article{Andreas:2012mt,
 author         = "Andreas, Sarah and Niebuhr, Carsten and Ringwald, Andreas",
 title          = "{New Limits on Hidden Photons from Past Electron 
                   Beam Dumps}",
 eprint         = "1209.6083",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "DESY-12-054",
 doi            = "10.1103/PhysRevD.86.095019",
 journal        = "Phys. Rev. D",
 volume         = "86",
 pages          = "095019",
 year           = "2012"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("limits/E137_Andreas2012mt.lmt")
efficiency = darkcast.Efficiency(lratio = 204.0/179.0)
