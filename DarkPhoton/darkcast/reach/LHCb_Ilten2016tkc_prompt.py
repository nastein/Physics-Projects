# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """

This limit is a projection for LHCb searches using an inclusive
di-muon final. This limit was
extracted from figure 1 (blue lines) of Ilten:2016tck.

The limit is prompt, so t0 = 0 and the definition of the maximum
proper lifetime from Ilten:2016tck is used, 0.04/(m - 2m_mu) + 0.1 ps.

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
bounds = darkcast.Dataset("reach/LHCb_Ilten2016tkc_prompt.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: (0.004/(m - 2*darkcast.pars.mfs['mu']) + 0.1)*1e-12)
