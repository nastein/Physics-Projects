# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit cannot be used for recasting, but is provided for
reference. Note that this limit does not use the minimal dark photon
model and makes strong assumptions about the decay of the Higgs. Four
prompt bounds are associated with this limit:
'ATLAS_Aad2015sms_br5.lmt', 'ATLAS_Aad2015sms_br10.lmt',
'ATLAS_Aad2015sms_br20.lmt', 'ATLAS_Aad2015sms_br40.lmt'. Each
corresponds to an assumed branching of the Higgs into two dark
photons: 5%, 10%, 20%, and 40% respectively. These were extracted from
figure 13 (red, orange, yellow, and green fills labeled ATLAS prompt)
of Aad:2015sms.
"""
bibtex = """
@article{Aad:2015sms,
 author         = "Aad, Georges and others",
 title          = "{A search for prompt lepton-jets in $pp$ collisions at
                   $\sqrt{s}=$ 8 TeV with the ATLAS detector}",
 collaboration  = "ATLAS",
 journal        = "JHEP",
 volume         = "02",
 year           = "2016",
 pages          = "062",
 doi            = "10.1007/JHEP02(2016)062",
 eprint         = "1511.05542",
 archivePrefix  = "arXiv",
 primaryClass   = "hep-ex",
 reportNumber   = "CERN-PH-EP-2015-242",
 SLACcitation   = "%%CITATION = ARXIV:1511.05542;%%"
}
"""
