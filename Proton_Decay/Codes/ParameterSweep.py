import numpy as np
import matplotlib.pyplot as plt
import math
pi = np.pi

hbar = 6.582122e-25
piplus_mass = .13957
pi0_mass = .13497
eta_mass = .5479
proton_mass = .938272
neutron_mass = .939565
p_pinubar_limit = 3.9e32
p_pieplus_limit = 1.7e34
p_etaeplus_limit = 1e34
V = .9743


def Phase_space(nucleon_mass, hadron_mass):

	return ((86400*365/hbar) * nucleon_mass/(32*pi) * ((1 - (hadron_mass/nucleon_mass)**2)**2))

def getmaxCeff(nucleon_mass,hadron_mass,Lambda,lifetime):

	return ((Lambda**4)/(Phase_space(nucleon_mass,hadron_mass)*lifetime))**.5

def Wilson_constraint_piplus_nubar(Lambda,lifetime):

	#Get the maximum value of Ceff based on the p-> e+ pi0 lifetime and Lambda
	maxCeff = getmaxCeff(proton_mass,pi0_mass,Lambda,p_pieplus_limit)

	#Find C1, C3 that give the above maximum Ceff
	c1 = np.logspace(-4,4, 100)
	c3 = (maxCeff - .131*c1*1.247*V)/(V*.134*1.25)

	good_c1 = []
	good_c3 = []

	#Get values of C1 and C3 from abov which ALSO satisfy p->pi+ nubar constraint
	for x,y in zip(c1,c3):
		if np.log10(((Lambda)**4)/(Phase_space(proton_mass,piplus_mass)*((-x*1.247*V*(-.186) + .189*(V**2)*1.25*y)**2))) > np.log10(lifetime):
			good_c1.append(x)
			good_c3.append(y)

	maxC1 = np.amax(good_c1)
	
	return good_c1,good_c3,maxC1

def Wilson_constraint_eta_eplus(Lambda,lifetime):

	#Get the maximum value of Ceff based on the p-> e+ pi0 lifetime and Lambda
	maxCeff = getmaxCeff(proton_mass,pi0_mass,Lambda,p_pieplus_limit)

	#Find C1, C3 that give the above maximum Ceff
	c1 = np.logspace(-4,4, 100)
	c3 = (maxCeff - .131*c1*1.247*V)/(V*.134*1.25)

	good_c1 = []
	good_c3 = []

	#Get values of C1 and C3 from abov which ALSO satisfy p->pi+ nubar constraint
	for x,y in zip(c1,c3):
		if np.log10(((Lambda)**4)/(Phase_space(proton_mass,eta_mass)*((-x*1.247*V*(.006) + .113*(V**2)*1.25*y)**2))) > np.log10(lifetime):
			good_c1.append(x)
			good_c3.append(y)

	maxC1 = np.amax(good_c1)
	
	return good_c1,good_c3,maxC1


"""
C1_16,C3_16,max_16 = Wilson_constraint_piplus_nubar(1e16,p_pinubar_limit)
C1_15,C3_15,max_15 = Wilson_constraint_piplus_nubar(1e15,p_pinubar_limit)
C1_14,C3_14,max_14 = Wilson_constraint_piplus_nubar(1e14,p_pinubar_limit)
C1_13,C3_13,max_13 = Wilson_constraint_piplus_nubar(1e13,p_pinubar_limit)

etaC1_16,etaC3_16,etamax_16 = Wilson_constraint_eta_eplus(1e16,p_etaeplus_limit)
etaC1_15,etaC3_15,etamax_15 = Wilson_constraint_eta_eplus(1e15,p_etaeplus_limit)
etaC1_14,etaC3_14,etamax_14 = Wilson_constraint_eta_eplus(1e14,p_etaeplus_limit)
#C1_13,C3_13,max_13 = Wilson_constraint_eta_eplus(1e13,p_etaeplus_limit)
"""


"""
lambdas = [1e14,1e15,1e16]
cmaxs = [max_14, max_15, max_16]
etacmaxs = [etamax_14, etamax_15, etamax_16]
plt.scatter(lambdas,cmaxs, c='red', label='p->pi+ nubar channel')
plt.scatter(lambdas, etacmaxs, c='blue', label='label = p->e+ eta channel')
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-5,1e4])
plt.xlim([5e11,2e16])
plt.xlabel('Lambda (GeV)')
plt.ylabel('Maximum value of C1 consistent with p->e+ pi0 channel')
plt.legend()
plt.show()
"""

"""
sc = plt.scatter(C1_14,C3_14)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([5e-5, 2])
#plt.ylim([5e-5, 2])
plt.xlabel('C1 values which allow for p->e+ pi0 cancellation')
plt.ylabel('-C3 values which allow for p->e+ pi0 cancellation')
plt.title('Log(p->pi+ v lifetime) for Lambda BNV = 1e14 GeV, satisfying p->e+ pi0 constraints')
plt.show()
"""

def Wilson_constraint_WET(Lambda, lifetime):

	#Get the maximum value of Ceff based on the p-> e+ pi0 lifetime and Lambda
	maxCeff = getmaxCeff(proton_mass,pi0_mass,Lambda,p_pieplus_limit)

	#Find C1, C3 that give the above maximum Ceff
	c1 = np.logspace(-7,3, 200)
	c3 = (maxCeff + .131*c1)/(.134)

	good_c1 = []
	good_c3 = []

	#Get values of C1 and C3 from abov which ALSO satisfy p->pi+ nubar constraint
	for x,y in zip(c1,c3):
		if np.log10(((Lambda)**4)/(Phase_space(proton_mass,piplus_mass)*((x*(.006) + .113*y)**2))) > np.log10(lifetime):
			good_c1.append(x)
			good_c3.append(y)

	maxC1 = np.amax(good_c1)
	maxC3 = (maxCeff + .131*maxC1)/(.134)
	
	return good_c1,good_c3,maxC1,maxC3

C1_16,C3_16,maxC1_16,maxC3_16 = Wilson_constraint_WET(1e16,p_etaeplus_limit)
C1_15,C3_15,maxC1_15,maxC3_15 = Wilson_constraint_WET(1e15,p_etaeplus_limit)
C1_14,C3_14,maxC1_14,maxC3_14 = Wilson_constraint_WET(1e14,p_etaeplus_limit)
C1_13,C3_13,maxC1_13,maxC3_13 = Wilson_constraint_WET(1e13,p_etaeplus_limit)

lambdas = [1e13,1e14,1e15,1e16]
c1maxs = [maxC1_13, maxC1_14, maxC1_15, maxC1_16]
c3maxs = [maxC3_13, maxC3_14, maxC3_15, maxC3_16]

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.scatter(lambdas, c1maxs, c='blue')
ax2.scatter(lambdas, c3maxs, c='red')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax1.set_ylim([5e-8,5e3])
ax2.set_ylim([5e-8,5e3])
ax1.set_xlim([5e12,2e16])
ax1.set_xlabel('Lambda (GeV)')
ax1.set_ylabel('Maximum value of C1 consistent with p->e+ pi0 lower limit', color='blue')
ax2.set_ylabel('Maximum value of C3 consistent with p->e+ pi0 lower limit', color='red')
ax1.set_title('Wilson Coefficients C1 and C3 vs BNV Scale for p->e+ eta channel')
plt.show()













