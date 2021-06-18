# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit is a projection for HPS searches using a di-electron final
state and was provided by the authors.

The limit is displaced, but does not have r-values and is not a beam
dump, and so the extrapolation behaviour for displaced r-values is
used.
"""
bibtex = """
@article{Baltzell:2016eee,
 author         = "Baltzell, N. and others",
 title          = "{The Heavy Photon Search beamline and its performance}",
 collaboration  = "HPS",
 journal        = "Nucl. Instrum. Meth.",
 volume         = "A859",
 year           = "2017",
 pages          = "69-75",
 doi            = "10.1016/j.nima.2017.03.061",
 eprint         = "1612.07821",
 archivePrefix  = "arXiv",
 primaryClass   = "physics.ins-det",
 reportNumber   = "JLAB-PHY-17-2393",
 SLACcitation   = "%%CITATION = ARXIV:1612.07821;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production("e_brem")
decay = "e_e"
bounds = darkcast.Datasets("reach/HPS_Baltzell2016eee_displaced.lmt")
efficiency = darkcast.Efficiency(lratio = float("inf"))
