# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """

This limit is a projection for LHCb searches using an inclusive
di-muon final. This limit was
extracted from figure 1 (blue lines) of Ilten:2016tck.

The limit is displaced, but does not have r-values and is not a beam
dump, and so the extrapolation behaviour for displaced r-values is
used.

The production is non-trivial and so the Monte Carlo results of
Ilten:2018crw are used, stored in 'LHCb_Aaij2017rft.prd'.
"""
bibtex = """
@article{Ilten:2016tkc,
 author         = "Ilten, Philip and Soreq, Yotam and Thaler, Jesse and
                   Williams, Mike and Xue, Wei",
 title          = "{Proposed Inclusive Dark Photon Search at LHCb}",
 journal        = "Phys. Rev. Lett.",
 volume         = "116",
 year           = "2016",
 number         = "25",
 pages          = "251803",
 doi            = "10.1103/PhysRevLett.116.251803",
 eprint         = "1603.08926",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ph",
 reportNumber   = "MIT-CTP-4785",
 SLACcitation   = "%%CITATION = ARXIV:1603.08926;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production(
    darkcast.Datasets("limits/LHCb_Aaij2017rft.prd"))
decay = "mu_mu"
bounds = darkcast.Datasets("reach/LHCb_Ilten2016tkc_pre.lmt")
efficiency = darkcast.Efficiency(lratio = float("inf"))
