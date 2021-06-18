# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit was taken from figure 3 of CMS:2019kiy.

The limit is prompt, so t0 = 0. For dark photon masses below 45 GeV
there are no lifetime requirements and so t1 = infinity. Above 45 GeV
transverse and longitudinal requirements are placed on the impact
parameters of the muons. The average maximum lifetime for these
requirements is approximated from Pythia 8 simulation and found to be
rougly 5e-12 seconds.

The non-trivial production is assumed to be primarily from Drell-Yan
production, is simulated with Pythia 8, and is stored in
'CMS_CMS2019kiy.prd'.
"""
bibtex = """
@article{CMS:2019kiy,
 author         = "CMS Collaboration",
 title          = "{Search for a narrow resonance decaying to a pair of
                   muons in proton-proton collisions at 13 TeV}",
 collaboration  = "CMS",
 year           = "2019",
 reportNumber   = "CMS-PAS-EXO-19-018",
 SLACcitation   = "%%CITATION = CMS-PAS-EXO-19-018;%%"
}
"""
model = darkcast.Model("dark_photon")
production = darkcast.Production(
    darkcast.Datasets("limits/CMS_CMS2019kiy.prd"))
decay = "mu_mu"
bounds = darkcast.Dataset("limits/CMS_CMS2019kiy.lmt")
def t1(m):
    if m < 45:  return float("inf")     # No IP requirements.
    if m < 130: return 5e-12            # Roughly flat lifetime.
    else: return -3.66e-11 + 3.12e-13*m # Roughly linear lifetime.
efficiency = darkcast.Efficiency(t0 = 0, t1 = t1)
