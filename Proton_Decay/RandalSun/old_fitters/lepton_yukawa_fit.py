import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.optimize
from scipy.optimize import fsolve
import sys
import csv
from numpy import sqrt as sqrt
from numpy import log as log
pi = np.pi

#Kr as a function of Lambda_IR
def kr(lam):
	return np.log((2.435e18)/lam)/pi

#quark yukawa coupling as a function of cl,cr, and Lambda_IR
def qyuk(cl,cr,lam):
	return sqrt((np.exp(2*(1 - cl - cr)*pi*kr(lam))*(.5 - cl)*(.5 - cr))/((np.exp((1 - 2*cl)*pi*kr(lam)) - 1)*(np.exp((1 - 2*cr)*pi*kr(lam)) - 1)))

me = []
mmu = []
mtau = []
lambdas = np.logspace(3,17,60)


j=0
with open('lepton_masses.csv','r') as h:
	reader = csv.reader(h)
	for row in reader:
		if j == 0:
			me = np.asfarray(np.array(row),float)
		if j == 1:
			mmu = np.asfarray(np.array(row),float)
		if j == 2:
			mtau = np.asfarray(np.array(row),float)
		j = j+1

yu = []
yd = []
ye = []


for i in range(0,len(lambdas)):
	ye.append([[me[i], 0., 0.],[0., mmu[i], 0.],[0., 0., mtau[i]]])

def calc_c_e(lam,ye):
	def f(c):
		return np.asarray((
			qyuk(c,c,lam) - ye[0][0]
			))
	
	#Construct an initial guess
	c0 = np.empty(1)
	c0.fill(.7)

	#Get the solution for the system
	c = scipy.optimize.least_squares(f,c0,method='lm').x
	return c[0]

def calc_c_mu(lam,ye):
	def f(c):
		return np.asarray((
			qyuk(c,c,lam) - ye[1][1]
			))
	
	#Construct an initial guess
	c0 = np.empty(1)
	c0.fill(.7)

	#Get the solution for the system
	c = scipy.optimize.least_squares(f,c0,method='lm').x
	return c[0]

def calc_c_tau(lam,ye):
	def f(c):
		return np.asarray((
			qyuk(c,c,lam) - ye[2][2]
			))
	
	#Construct an initial guess
	c0 = np.empty(1)
	c0.fill(.7)

	#Get the solution for the system
	c = scipy.optimize.least_squares(f,c0,method='lm').x
	return c[0]

elec_c = []
mu_c = []
tau_c = []
pred_electron_yukawa = []
pred_muon_yukawa = []
pred_tau_yukawa = []
electron_res = []
muon_res = []
tau_res = []

for i in range(0,len(lambdas)):
	ce = calc_c_e(lambdas[i],ye[i])
	cmu = calc_c_mu(lambdas[i],ye[i])
	ctau = calc_c_tau(lambdas[i],ye[i])
	elec_c.append(ce)
	mu_c.append(cmu)
	tau_c.append(ctau)
	pred_electron_yukawa.append(qyuk(ce,ce,lambdas[i]))
	pred_muon_yukawa.append(qyuk(cmu,cmu,lambdas[i]))
	pred_tau_yukawa.append(qyuk(ctau,ctau,lambdas[i]))
	electron_res.append(qyuk(ce,ce,lambdas[i]) - me[i])
	muon_res.append(qyuk(cmu,cmu,lambdas[i]) - mmu[i])
	tau_res.append(qyuk(ctau,ctau,lambdas[i]) - mtau[i])

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.set_xscale('log')
ax1.scatter(lambdas,elec_c,c='red',label='electron c')
ax1.scatter(lambdas,mu_c,c='brown',label='muon c')
ax1.scatter(lambdas,tau_c,c='black',label='tau c')
ax1.legend(loc=3)
ax1.set_title('5d mass parameters for left handed charged leptons')

ax2 =fig.add_subplot(212)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.scatter(lambdas,pred_electron_yukawa,c='red',label='electron yukawa')
ax2.scatter(lambdas,pred_muon_yukawa,c='brown',label='muon yukawa coupling')
ax2.scatter(lambdas,pred_tau_yukawa,c='black',label='tau yukawa coupling')
ax2.legend(loc=3)
ax2.set_title('Predicted 4d effective yukawa coupling')

np.save('elec_c',elec_c)
np.save('muon_c',mu_c)
np.save('tau_c',tau_c)

plt.show()









