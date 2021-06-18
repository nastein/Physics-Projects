# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was taken from the 'long-lived.txt' data file supplied in the
supplemental material of LHCb_Aaij2017rft, available at:
https://cds.cern.ch/record/2287638/files/LHCb-PAPER-2017-038-figures.zip

The limit is displaced, with the the upper limit on the observed dark
photon yield relative to the expected number of dark photon decays
given for each mass and global coupling. Outside the provided limits
the efficiency is approximated with t0 = 0 and t1 = tau_max for
lifetimes above the maximum, while t0 = tau_min and t1 = infinity for
lifteimes below the minimum.

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
bounds = darkcast.Dataset("limits/LHCb_Aaij2017rft_displaced.lmt")
efficiency = darkcast.Efficiency(rvals = True)
