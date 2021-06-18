# DARKCAST is licensed under the GNU GPL version 2 or later.
# Copyright (C) 2020 Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.
import darkcast
notes = """
This limit cannot be used for recasting, but is provided for
reference. Note that this limit does not use the minimal dark photon
model and makes strong assumptions about the decay of the Higgs. Four
prompt bounds are associated with this limit:
'CMS_CMS2016tgd_br1.lmt', 'CMS_CMS2016tgd_br5.lmt',
'CMS_CMS2016tgd_br10.lmt', 'CMS_CMS2016tgd_br40.lmt'. Each corresponds
to an assumed branching of the Higgs into two dark photons: 1%, 5%,
10%, and 40% respectively. These were extracted from figure 6 (shaded
orange fills labeled BR = 1, 5, 10, and 40%) of CMS:2016tgd.
"""
bibtex = """
@article{CMS:2016tgd,
 author         = "CMS Collaboration",
 title          = "{A Search for Beyond Standard Model Light Bosons Decaying
                   into Muon Pairs}",
 collaboration  = "CMS",
 year           = "2016",
 reportNumber   = "CMS-PAS-HIG-16-035",
 SLACcitation   = "%%CITATION = CMS-PAS-HIG-16-035;%%"
}
"""
