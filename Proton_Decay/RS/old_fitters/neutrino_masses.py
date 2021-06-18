import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.optimize
from scipy.optimize import fsolve
import sys
import csv
from numpy import sqrt as sqrt
from numpy import log as log
from numpy import exp as exp
pi = np.pi

m_pl = 2.435e18
v=246

#Kr as a function of Lambda_IR
def kr(lam):
	return np.log((m_pl)/lam)/pi

elec_c = np.load('elec_c.npy')
muon_c = np.load('muon_c.npy')
tau_c = np.load('tau_c.npy')
elec_k_c = np.load('elec_c_k.npy')
muon_k_c = np.load('muon_c_k.npy')
tau_k_c = np.load('tau_c_k.npy')

lambdas = np.logspace(3,17,60)

def W_scale(c,lam,k):
	return (k/(m_pl**2))*(.5 - c)*((exp(2*pi*kr(lam)*(1 - c)))/(exp(pi*kr(lam)*(1 - 2*c)) - 1))

W_elec_array = []
W_muon_array = []
W_tau_array = []
W_elec_k_array = []
W_muon_k_array = []
W_tau_k_array = []

for i,lam in enumerate(lambdas):
	W_elec_array.append(W_scale(elec_c[i],lam,m_pl)*(v**2)*1e9)
	W_muon_array.append(W_scale(muon_c[i],lam,m_pl)*(v**2)*1e9)
	W_tau_array.append(W_scale(tau_c[i],lam,m_pl)*(v**2)*1e9)
	W_elec_k_array.append(W_scale(elec_k_c[i],lam,.1*m_pl)*(v**2)*1e9)
	W_muon_k_array.append(W_scale(muon_k_c[i],lam,.1*m_pl)*(v**2)*1e9)
	W_tau_k_array.append(W_scale(tau_k_c[i],lam,.1*m_pl)*(v**2)*1e9)


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Lambda_IR (GeV)')
ax1.set_ylabel('neutrino mass (eV)')
ax1.scatter(lambdas,W_elec_array,c='red',label=r'$v^2(L_e L_e)/\Lambda_W$, $k=M_{pl}$')
ax1.scatter(lambdas,W_elec_k_array,c='purple',label=r'$v^2(L_e L_e)/\Lambda_W$, $k=\frac{M_{pl}}{10}$')
ax1.scatter(lambdas,W_muon_array,c='blue',label=r'$v^2(L_{\mu} L_{\mu})/\Lambda_W$, $k=M_{pl}$')
ax1.scatter(lambdas,W_muon_k_array,c='brown',label=r'$v^2(L_{\mu} L_{\mu})/\Lambda_W$, $k=\frac{M_{pl}}{10}$')
ax1.scatter(lambdas,W_tau_array,c='black',label=r'$v^2(L_{\tau}L_{\tau})/\Lambda_W$, $k=M_{pl}$')
ax1.scatter(lambdas,W_tau_k_array,c='orange',label=r'$v^2(L_{\tau}L_{\tau})/\Lambda_W$, $k=\frac{M_{pl}}{10}$')
ax1.legend(loc=3)
ax1.set_title('Neutrino masses from D=5 Weinberg operator')

plt.show()





