import numpy as np
import pylab as py
import matplotlib.pyplot as plt 
import math
import time
import random as rand

#Simulate the ising model
def compute_energy(lat, J):
	E = -J*(np.sum(np.roll(lat, 1, 0)*lat) + np.sum(np.roll(lat, 1, 1)*lat))
	return E

def compute_magnetization(lat,N):
	M = np.sum(lat)/N*N
	return np.abs(M)


def compute_site_energy(lat,J,x,y):
	E = -J*(lat[x,y]*lat[x-1,y] +  lat[x,y]*lat[np.mod(x+1,N),y] + lat[x,y]*lat[x,np.mod(y+1,N)] + lat[x,y]*lat[x,y-1])
	return E

def spin_flip(lat,N,J,T):
	x=rand.randint(0,N-1)
	y=rand.randint(0,N-1)

	new_lat = np.copy(lat)
	E_orig = compute_site_energy(lat, J, x, y)
	
	new_lat[x,y] = lat[x,y]*(-1)
	E_new = compute_site_energy(new_lat, J, x, y)
	E_diff = E_new - E_orig

	if (E_diff < 0):
		return 1, new_lat, E_diff
	else:
		r = np.random.random()
		if r < np.exp(-E_diff*(1/T)):
			return 1, new_lat, E_diff
		else:
			return 0, lat, E_diff


#Create lattice and initialize coupling constant
J = 1
N = 10

lat = np.ones((N,N))
for i in range(0,N):
	for j in range(0,N):
		u = np.random.random()
		if u < .5:
			lat[i,j] = -1

#Create array of temperatures and average observables
temps = np.arange(1, 11, .01)
Energy = np.zeros(len(temps))
EnergySq = np.zeros(len(temps))
Magnetization = np.zeros(len(temps))
MagnetizationSq = np.zeros(len(temps))
suscepitability = np.zeros(len(temps))
heat_capacity = np.zeros(len(temps))
Z = np.zeros(len(temps))
lattice = np.copy(lat)

original_E = compute_energy(lat,J)
E=original_E

for t in range(0,len(temps)):

	for i in range (0, 6000):
		if (i == 0):
			Energy[t] += compute_energy(lat,J)*np.exp(-original_E*(1/temps[t]))
			EnergySq[t] +=(compute_energy(lat,J)**2)*np.exp(-original_E*(1/temps[t]))
			Magnetization[t] += compute_magnetization(lat,N)*np.exp(-original_E*(1/temps[t]))
			MagnetizationSq[t] += (compute_magnetization(lat,N)**2)*np.exp(-original_E*(1/temps[t]))
			Z[t] += np.exp(-original_E*(1/temps[t]))
		else:
			flipped, lattice, Edif = spin_flip(lattice,N,J,temps[t])
			if (flipped == 0):
				continue
			else:
				E = E+Edif
				Energy[t] += E*np.exp(-E*(1/temps[t]))
				EnergySq[t] +=(E**2)*np.exp(-E*(1/temps[t]))
				Magnetization[t] += compute_magnetization(lattice,N)*np.exp(-E*(1/temps[t]))
				MagnetizationSq[t] += (compute_magnetization(lattice,N)**2)*np.exp(-E*(1/temps[t]))
				Z[t] += np.exp(-E*(1/temps[t]))

	Energy[t] = Energy[t]/Z[t]
	EnergySq[t] = EnergySq[t]/Z[t]
	Magnetization[t] = Magnetization[t]/Z[t]
	MagnetizationSq[t] = MagnetizationSq[t]/Z[t]

for i in range (0,len(temps)):
	heat_capacity[i] = (EnergySq[i] - (Energy[i]**2))/(temps[i]**2)
	suscepitability[i] = (MagnetizationSq[i] - (Magnetization[i]**2))/temps[i]



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(temps, Energy)
ax1.set_xlabel('Temperature (k_b = 1)')
ax1.set_ylabel('Average Energy ')
ax1.set_title('Average Energy vs Temperature, J = 1')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(temps, Magnetization)
ax2.set_xlabel('Temperature (k_b = 1)')
ax2.set_ylabel('Average Magnetization')
ax2.set_title('Average Magnetization vs Temperature, J = 1')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(temps, suscepitability)
ax3.set_xlabel('Temperature (k_b = 1)')
ax3.set_ylabel('Average Suscepitability')
ax3.set_title('Average Suscepitability vs Temperature, J = 1')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(temps, heat_capacity)
ax4.set_xlabel('Temperature (k_b = 1)')
ax4.set_ylabel('Average Heat Capacity')
ax4.set_title('Average Heat Capacity vs Temperature, J = 1')

plt.show()






