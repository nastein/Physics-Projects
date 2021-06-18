# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit cannot be used for recasting, but is provided for
reference. Note that this limit does not use the minimal dark photon
model and makes strong assumptions about the decay of the Higgs. Six
diplaced bounds are associated with this limit:
'ATLAS_Aad2014yea_br5[ab].lmt', 'ATLAS_Aad2014yea_br10[ab].lmt',
'ATLAS_Aad2014yea_br20.lmt', 'ATLAS_Aad2014yea_br40.lmt'. Each
corresponds to an assumed branching of the Higgs into two dark
photons: 5%, 10%, 20%, and 40% respectively. When an 'a' and 'b' bound
are given these constitute two independent exclusion regions. All
bounds were extracted from figure 13 (red, orange, yellow, and green
fills labeled ATLAS displaced) of Aad:2015sms.
"""
bibtex = """
@article{Aad:2014yea,
 author         = "Aad, Georges and others",
 title          = "{Search for long-lived neutral particles decaying into
                   lepton jets in proton-proton collisions at $ \sqrt{s}=8 $
                   TeV with the ATLAS detector}",
 collaboration  = "ATLAS",
 journal        = "JHEP",
 volume         = "11",
 year           = "2014",
 pages          = "088",
 doi            = "10.1007/JHEP11(2014)088",
 eprint         = "1409.0746",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-PH-EP-2014-209",
 SLACcitation   = "%%CITATION = ARXIV:1409.0746;%%"
}
"""
