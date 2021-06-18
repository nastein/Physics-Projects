import numpy as np
import pylab as py
import matplotlib.pyplot as plt 
import math
import time
import random as rand

############ 3D Ising model simulation with Swednson Wang clustering ################

def init_lattice(N):
	lat = np.ones((N,N,N))
	for i in range (0,N):
		for j in range (0,N):
			for k in range (0,N):
				u = np.random.random()
				if (u <= .5):
					lat[i,j,k] = -1
				else:
					lat[i,j,k] = 1
	return lat

def compute_energy(lat, J,N):
	E = -(J)*(np.sum(np.roll(lat, 1, 0)*lat) + np.sum(np.roll(lat, 1, 1)*lat) + np.sum(np.roll(lat, 1, 2)*lat))
	return (E)

def compute_magnetization(lat,N):
	M = np.sum(lat)
	return np.abs(M)

def freeze_probability(J,t):
	return 1 - np.exp(-2*J/(t))

def findEC(n,proper_labels):
	m = n
	while proper_labels[m] != m:
		m = proper_labels[m]
	while proper_labels[n] != n:
		l = proper_labels[n]
		proper_labels[n] = m
		n = l
	return m

def unionEC(m,n,proper_labels):
	proper_labels[findEC(m,proper_labels)] = findEC(n,proper_labels)

def relabel(proper_labels,lbarray,N):
	n = 0
	for i in np.arange(0,N):
		for j in np.arange(0,N):
			for k in np.arange(0,N):
				n =  n + 1
				if lbarray[i,j,k] != 0:
					y = lbarray[i,j,k]
					while proper_labels[y] != y:
						y = proper_labels[y]
					lbarray[i,j,k] = y
	return lbarray

def make_bonds(N,lat, J, t):
	labels = np.zeros([N,N,N], dtype = int)
	prop_labels = [0]
	latest_label = 0

	freezeProbability = freeze_probability(J, t)

	for i in range (0,N):
		for j in range (0,N):
			for k in range (0,N):
				
				current_spin = lat[i,j,k]

				if (i!=0):
					above = lat[i-1,j,k]
				else:
					above = 0
				if (j!=0):
					left = lat[i,j-1,k]
				else:
					left = 0
				if (k!=0):
					behind = lat[i,j,k-1]
				else:
					behind = 0

				if (above != current_spin and (left != current_spin) and (behind != current_spin)):
					latest_label = latest_label + 1
					labels[i,j,k] = latest_label
					prop_labels.append(latest_label)

				elif (above == current_spin and (left != current_spin) and (behind != current_spin)):
					u = np.random.random()
					if (u <= freezeProbability):
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)
					else:
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)

				elif (above != current_spin and (left == current_spin) and (behind != current_spin)):
					u = np.random.random()
					if (u <= freezeProbability):
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)
					else:
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)

				elif (above != current_spin and (left != current_spin) and (behind == current_spin)):
					u = np.random.random()
					if (u <= freezeProbability):
						labels[i,j,k] = findEC(labels[i,j,k-1], prop_labels)
					else:
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)
						

				elif (above == current_spin and (left == current_spin) and (behind != current_spin)):
					u = np.random.random()
					v = np.random.random()

					if (u <= freezeProbability and v > freezeProbability):
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					elif (u > freezeProbability and v <= freezeProbability):
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)

					elif (u <= freezeProbability and v <= freezeProbability):	
						unionEC(labels[i-1,j,k], labels[i,j-1,k], prop_labels)
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					elif(u > freezeProbability and v > freezeProbability):
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)
						

				elif (above == current_spin and (left != current_spin) and (behind == current_spin)):
					u = np.random.random()
					v = np.random.random()

					if (u <= freezeProbability and v > freezeProbability):
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					elif (u > freezeProbability and v <= freezeProbability):
						labels[i,j,k] = findEC(labels[i,j,k-1], prop_labels)

					elif (u <= freezeProbability and v <= freezeProbability):	
						unionEC(labels[i-1,j,k], labels[i,j,k-1], prop_labels)
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					elif(u > freezeProbability and v > freezeProbability):
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)

				elif (above != current_spin and (left == current_spin) and (behind == current_spin)):
					u = np.random.random()
					v = np.random.random()

					if (u <= freezeProbability and v > freezeProbability):
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)

					elif (u > freezeProbability and v <= freezeProbability):
						labels[i,j,k] = findEC(labels[i,j,k-1], prop_labels)

					elif (u <= freezeProbability and v <= freezeProbability):	
						unionEC(labels[i,j-1,k], labels[i,j,k-1], prop_labels)
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)

					elif(u > freezeProbability and v > freezeProbability):
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)


				elif (above == current_spin and (left == current_spin) and (behind == current_spin)):
					u = np.random.random()
					v = np.random.random()
					w = np.random.random()

					if (u > freezeProbability and v > freezeProbability and w > freezeProbability):
						latest_label = latest_label + 1
						labels[i,j,k] = latest_label
						prop_labels.append(latest_label)

					if (u <= freezeProbability and v > freezeProbability and w > freezeProbability):
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					if (u > freezeProbability and v <= freezeProbability and w > freezeProbability):
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)

					if (u > freezeProbability and v > freezeProbability and w <= freezeProbability):	
						labels[i,j,k] = findEC(labels[i,j,k-1], prop_labels)

					if (u <= freezeProbability and v <= freezeProbability and w > freezeProbability):	
						unionEC(labels[i-1,j,k], labels[i,j-1,k], prop_labels)
						labels[i,j,k] = findEC(labels[i,j,k-1], prop_labels)

					if (u <= freezeProbability and v > freezeProbability and w <= freezeProbability):	
						unionEC(labels[i-1,j,k], labels[i,j,k-1], prop_labels)
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)

					if (u > freezeProbability and v <=freezeProbability and w <= freezeProbability):	
						unionEC(labels[i,j-1,k], labels[i,j,k-1], prop_labels)
						labels[i,j,k] = findEC(labels[i,j-1,k], prop_labels)

					if (u <= freezeProbability and v <= freezeProbability and w <= freezeProbability):	
						unionEC(labels[i-1,j,k], labels[i,j,k-1], prop_labels)
						unionEC(labels[i,j-1,k], labels[i,j,k-1], prop_labels)
						unionEC(labels[i-1,j,k], labels[i,j-1,k], prop_labels)
						labels[i,j,k] = findEC(labels[i-1,j,k], prop_labels)


	labels = relabel(prop_labels, labels, N)

	for i in range(N):
		for j in range(N):
			pboundary = np.random.rand()
			if pboundary <= freezeProbability and lat[i,j,0] == lat[i,j,-1]:
				unionEC(labels[i,j,0],labels[i,j,-1],prop_labels)
				#label_array[i,j,k] = findEC(label_array[i,j-1,k],proper_labels)
	for i in range(N):
		for k in range(N):
			pboundary = np.random.rand()
			if pboundary <= freezeProbability and lat[i,0,k] == lat[i,-1,k]:
				unionEC(labels[i,0,k],labels[i,-1,k],prop_labels)
				#label_array[i,j,k] = findEC(label_array[i,j-1,k],proper_labels)
	for j in range(N):
		for k in range(N):
			pboundary = np.random.rand()
			if pboundary <= freezeProbability and lat[0,j,k] == lat[-1,j,k]:
				unionEC(labels[0,j,k],labels[-1,j,k],prop_labels)
				#label_array[i,j,k] = findEC(label_array[i,j-1,k],proper_labels)

	labels = relabel(prop_labels, labels, N)

	return labels, prop_labels


def flip_spins(lat,label_array,proper_labels,N):
	for l in range(0,proper_labels.size):
		u = np.random.random()
		indices = np.where(label_array == proper_labels[l])

		if u <= .5:
			lat[indices] = lat[indices]*(-1)
		else:
			continue

	return lat

def compute_observables(N, J, temps, therm_steps, mc_steps):

	print('Lattice size = ', N)

	#Create array of temperatures and average observables
	Energy = np.zeros(len(temps))
	EnergySq = np.zeros(len(temps))
	Magnetization = np.zeros(len(temps))
	MagnetizationSq = np.zeros(len(temps))
	MagnetizationFour = np.zeros(len(temps))
	susceptibility = np.zeros(len(temps))
	heat_capacity = np.zeros(len(temps))
	Z = np.zeros(len(temps))
	binder = np.zeros(len(temps))

	#lat = init_lattice(N)

	for t in range(0,len(temps)):
		print('Temperature = ', temps[t])

		#re initialize the lattice 
		lat = init_lattice(N)

		for i in range (0, mc_steps):

			if (i <= therm_steps):
				labels, proper_labels = make_bonds(N, lat, J, temps[t])
				final_proper_labels = np.unique(proper_labels)
				lat = flip_spins(lat, labels, final_proper_labels,N)

			else:
				labels, proper_labels = make_bonds(N, lat, J, temps[t])
				final_proper_labels = np.unique(proper_labels)
				lat = flip_spins(lat, labels, final_proper_labels,N)
				E = compute_energy(lat,J,N)

				Energy[t] += E
				EnergySq[t] += (E**2)
				Magnetization[t] += compute_magnetization(lat,N)
				MagnetizationSq[t] += compute_magnetization(lat,N)**2
				MagnetizationFour[t] += compute_magnetization(lat,N)**4
				Z[t] += 1

		Energy[t] = Energy[t]/Z[t]
		EnergySq[t] = EnergySq[t]/Z[t]
		Magnetization[t] = Magnetization[t]/Z[t]
		MagnetizationSq[t] = MagnetizationSq[t]/Z[t]
		MagnetizationFour[t] = MagnetizationFour[t]/Z[t]

	for i in range (0,len(temps)):
		heat_capacity[i] = (EnergySq[i] - (Energy[i]**2))/((N**3)*temps[i]**2)
		susceptibility[i] = ((MagnetizationSq[i]) - (Magnetization[i]**2))/(temps[i]*(N**3))
		binder[i] = 1 - ((MagnetizationFour[i])/(3*(MagnetizationSq[i]**2)))

	#Temperature where the susceptability is maximal (should be T_c)
	crit_temp_index_sus = np.squeeze(np.where(np.max(susceptibility) == susceptibility))
	crit_temp_index_heat = np.squeeze(np.where(np.max(heat_capacity) == heat_capacity))

	return Energy, Magnetization, heat_capacity, susceptibility, binder, crit_temp_index_sus, crit_temp_index_heat

###################################################################################################
lattice_size_array = [6,8,10,12]
J = +1
log_lattice_size_array = np.log(lattice_size_array)
Energies = []
Magnetizations = []
Heat_Capacities = []
Susceptabilities = []
Binders = []

log_peak_magnetization = []
log_peak_heat_capacity = []
log_peak_susceptability = []

t_0 = 4.4
t_f = 5
t_step_size = .2
#temps = np.arange(t_0, t_f, t_step_size)
temps = [4.35,4.40,4.45,4.50]
therm_steps = 1000
mc_steps = 5000

for size in lattice_size_array:

	energy, magnetization, heat_capacity, susceptibility, binder, crit_temp_index_sus, crit_temp_index_heat = compute_observables(size, J, temps, therm_steps, mc_steps)

	Energies.append(energy/(size**3))
	Magnetizations.append(magnetization/(size**3))
	Heat_Capacities.append(heat_capacity)
	Susceptabilities.append(susceptibility)
	Binders.append(binder)

	print('peak temperature derived from susceptibility for size N = ', size,' is T = ', temps[crit_temp_index_sus])
	print('peak temperature derived from heat capacity for size N = ', size, ' is T = ', temps[crit_temp_index_heat])

	log_peak_magnetization.append(np.log(magnetization[crit_temp_index_sus]/(size**3)))
	log_peak_heat_capacity.append(np.log(np.max(heat_capacity)))
	log_peak_susceptability.append(np.log(np.max(susceptibility)))


###################################################################################################

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for size in range(0,len(lattice_size_array)):
	ax1.plot(temps, Energies[size], 'o', label = "{}{}".format("N = ",lattice_size_array[size]))
ax1.set_xlabel('Temperature (k_b = 1)')
ax1.set_ylabel('Average Energy ')
ax1.set_title('Average Energy vs Temperature, J = 1')
legend = ax1.legend(loc='upper right', shadow=False, fontsize='medium')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for size in range(0,len(lattice_size_array)):
	ax2.plot(temps, Magnetizations[size], 'o', label = "{}{}".format("N = ",lattice_size_array[size]))
ax2.set_xlabel('Temperature (k_b = 1)')
ax2.set_ylabel('Average Magnetization')
ax2.set_title('Average Magnetization vs Temperature, J = 1')
legend = ax2.legend(loc='upper right', shadow=False, fontsize='medium')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
for size in range(0,len(lattice_size_array)):
	ax3.plot(temps, Susceptabilities[size], 'o', label = "{}{}".format("N = ",lattice_size_array[size]))
ax3.set_xlabel('Temperature (k_b = 1)')
ax3.set_ylabel('Average susceptibility')
ax3.set_title('Average susceptibility vs Temperature, J = 1')
legend = ax3.legend(loc='upper right', shadow=False, fontsize='medium')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for size in range(0,len(lattice_size_array)):
	ax4.plot(temps, Heat_Capacities[size], 'o', label = "{}{}".format("N = ",lattice_size_array[size]))
ax4.set_xlabel('Temperature (k_b = 1)')
ax4.set_ylabel('Average Heat Capacity')
ax4.set_title('Average Heat Capacity vs Temperature, J = 1')
legend = ax4.legend(loc='upper right', shadow=False, fontsize='medium')

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
for size in range(0,len(lattice_size_array)):
	ax5.plot(temps, Binders[size], label = "{}{}".format("N = ",lattice_size_array[size]))
ax5.set_xlabel('Temperature (k_b = 1)')
ax5.set_ylabel('Binders cumulant')
ax5.set_title('Binders cumulant vs Temperature, J = 1')
legend = ax5.legend(loc='upper right', shadow=False, fontsize='medium')

fit_mag = np.polyfit(log_lattice_size_array, log_peak_magnetization, 1)
mag_fit_fn = np.poly1d(fit_mag)
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(log_lattice_size_array, log_peak_magnetization, 'o', label='data')
ax6.plot(log_lattice_size_array, mag_fit_fn(log_lattice_size_array), label="{}{}".format("fit, slope = ", fit_mag[0]))
ax6.set_xlabel('ln(Lattice Size = N)')
ax6.set_ylabel('ln(magnetization)')
ax6.set_title('Finite scaling relation for magnetization')
legend = ax6.legend(loc='upper right', shadow=False, fontsize='medium')

fit_sus = np.polyfit(log_lattice_size_array, log_peak_susceptability, 1)
sus_fit_fn = np.poly1d(fit_sus)
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.plot(log_lattice_size_array, log_peak_susceptability, 'o', label='data')
ax7.plot(log_lattice_size_array, sus_fit_fn(log_lattice_size_array), label="{}{}".format("fit, slope = ", fit_sus[0]))
ax7.set_xlabel('ln(Lattice Size = N)')
ax7.set_ylabel('ln(susceptibility)')
ax7.set_title('Finite scaling relation for susceptibility')
legend = ax7.legend(loc='upper right', shadow=False, fontsize='medium')

fit_cap = np.polyfit(log_lattice_size_array, log_peak_heat_capacity, 1)
heat_fit_fn = np.poly1d(fit_cap)
fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.plot(log_lattice_size_array,log_peak_heat_capacity, 'o', label='data')
ax8.plot(log_lattice_size_array, heat_fit_fn(log_lattice_size_array), label="{}{}".format("fit, slope = ", fit_cap[0]))
ax8.set_xlabel('ln(Lattice Size = N)')
ax8.set_ylabel('ln(heat capacity)')
ax8.set_title('Finite scaling relation for heat capacity')
legend = ax8.legend(loc='upper right', shadow=False, fontsize='medium')


plt.show()
