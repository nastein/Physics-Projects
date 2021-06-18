# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was taken from the 'prompt-excluded.txt' data file supplied
in the supplemental material of LHCb_Aaij2019bvg, available at:
https://cds.cern.ch/record/2691573/files/LHCb-PAPER-2019-031-figures.zip

The limit is prompt, so t0 = 0. The selection is similar to the
phenomenology study of Ilten:2016tkc and so the definition of the
maximum proper lifetime from this paper is used, 0.04/(m - 2m_mu) +
0.1 ps.

The production is non-trivial and so the Monte Carlo results of
Ilten:2018crw are used, stored in 'LHCb_Aaij2017rft.prd'. Note that
differences in the fiducial requirements will introduce some changes
in these production fractions, so proceed with caution when recasting
these results.
"""
bibtex = """
@article{Aaij:2019bvg,
 author         = "Aaij, Roel and others",
 title          = "{Search for $A'\!\to\!\mu^+\mu^-$ decays}",
 collaboration  = "LHCb",
 year           = "2019",
 eprint         = "1910.06926",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "LHCb-PAPER-2019-031, CERN-EP-2019-212",
 SLACcitation   = "%%CITATION = ARXIV:1910.06926;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production(
    darkcast.Datasets("limits/LHCb_Aaij2017rft.prd"))
decay = "mu_mu"
bounds = darkcast.Dataset("limits/LHCb_Aaij2019bvg_prompt.lmt")
efficiency = darkcast.Efficiency(
    t0 = 0, t1 = lambda m: (0.004/(m - 2*darkcast.pars.mfs['mu']) + 0.1)*1e-12)
