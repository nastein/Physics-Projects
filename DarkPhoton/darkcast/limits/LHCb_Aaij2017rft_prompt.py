# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was taken from the 'prompt.txt' data file supplied in the
supplemental material of LHCb_Aaij2017rft, available at:
https://cds.cern.ch/record/2287638/files/LHCb-PAPER-2017-038-figures.zip

The limit is prompt, so t0 = 0. The selection is the same as the
phenomenology study of Ilten:2016tkc and so the definition of the
maximum proper lifetime from this paper is used, 0.04/(m - 2m_mu) +
0.1 ps.

The production is non-trivial and so the Monte Carlo results of
Ilten:2018crw are used, stored in 'LHCb_Aaij2017rft.prd'.
"""
bibtex = """
@article{Aaij:2017rft,
 author         = "Aaij, Roel and others",
 title          = "{Search for Dark Photons Produced in 13 TeV $pp$
                   Collisions}",
 collaboration  = "LHCb",
 journal        = "Phys. Rev. Lett.",
 volume         = "120",
 year           = "2018",
 number         = "6",
 pages          = "061801",
 doi            = "10.1103/PhysRevLett.120.061801",
 eprint         = "1710.02867",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "LHCB-PAPER-2017-038, CERN-EP-2017-248",
 SLACcitation   = "%%CITATION = ARXIV:1710.02867;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production(
    darkcast.Datasets("limits/LHCb_Aaij2017rft.prd"))
decay = "mu_mu"
bounds = darkcast.Dataset("limits/LHCb_Aaij2017rft_prompt.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: (0.004/(m - 2*darkcast.pars.mfs['mu']) + 0.1)*1e-12)
