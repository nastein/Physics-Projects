import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as sc 

#pvalue = 0.05 #for 95% exclusion
g_nom = 1e-3

def get_sig(bkgd, pvalue = 0.05):
	#Formula for p:
	# 	p = Q(b + 1, s+ b)
	added_bkgd = bkgd + 1
	
	#below function returns s + b 
	solution = np.array(sc.gammainccinv(added_bkgd, pvalue))

	#subtract off b to get s
	solution = solution - bkgd

	return solution;

def get_g_excl(S):
	masses = []
	if (S == 500):
		masses = [5 + 2*i for i in range(173)]
		bkgd_events = []
		with open("/Users/noah/Physics/James_Research/ILC_PhotonFusion/500GeV_LBL.txt") as bkgd_f:
			bkgd_events = [float(i) for i in bkgd_f]

		nominal_sig_events = np.array([])
		with open("/Users/noah/Physics/James_Research/ILC_PhotonFusion/500GeV_yield.txt") as nom_sig_f:
			nominal_sig_events = [float(i) for i in nom_sig_f]

		for i,n in enumerate(nominal_sig_events):
			if n == 0.0:
				bkgd_events.pop(i)
				nominal_sig_events.pop(i)
				masses.pop(i)

	if (S == 250):
		masses = [5 + 2*i for i in range(73)]
		bkgd_events = []
		with open("/Users/noah/Physics/James_Research/ILC_PhotonFusion/250GeV_LBL.txt") as bkgd_f:
			bkgd_events = [float(i) for i in bkgd_f]

		nominal_sig_events = np.array([])
		with open("/Users/noah/Physics/James_Research/ILC_PhotonFusion/250GeV_yield.txt") as nom_sig_f:
			nominal_sig_events = [float(i) for i in nom_sig_f]

		for i,n in enumerate(nominal_sig_events):
			if n == 0.0:
				bkgd_events.pop(i)
				nominal_sig_events.pop(i)
				masses.pop(i)
	
	needed_sig_events = np.array(get_sig(np.array(bkgd_events)))

	g_excl = np.sqrt(needed_sig_events/np.array(nominal_sig_events))*g_nom
	

	return g_excl,masses
