import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import sqrt as sqrt
from numpy import log as log
pi = np.pi

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

MGUT = 1e16
#mu_array = [14,15,16,17]
mu_array = [17]
colors = ['r','b','g','y','m']
markers = ['*', 'o']

Mpl = 2.435e18
c1 = -2
c2 = -6
c3 = 4

carray = [0,-.1,.1, -.5, .5]
#c = 1


fig = plt.figure()


def lambda_12(Mv,Mhc,M24t,mu,c):
  return -(14/3) + 112*log(Mv/mu) + 2*log(M24t/mu) - (8/6)*log(Mhc/mu) + c*(np.sqrt(8/25)*Mv/Mpl)*(c1-c2)*(48*np.pi*np.pi)

def lambda_23(Mv,Mhc,M24o,M24t,mu,c):
  return -2 + 21*log(Mv/mu) + (1/2)*log(Mhc*M24o/(mu*mu)) - 2*log(M24t/mu) + c*(np.sqrt(8/25)*Mv/Mpl)*(c2-c3)*(48*np.pi*np.pi)

sm12 = [293.51, 193.251, 92.977, -7.30987, -107.609, -207.919, -308.24, -408.571, -508.912]
sm23 = [351.668, 297.305, 243.049, 188.885, 134.802, 80.7879, 26.8361, -27.0608, -80.9081]
text = [r'$10^{10}$',r'$10^{11}$',r'$10^{12}$',r'$10^{13}$',r'$10^{14}$',r'$10^{15}$',r'$10^{16}$',r'$10^{17}$',r'$10^{18}$']

ax1 = fig.add_subplot(1,1,1)
ax1.scatter(sm12,sm23,color='black',zorder=2, label='SM threshold correction')

for i, txt in enumerate(text):
  ax1.annotate(txt, (sm12[i], sm23[i]))

for j,mu in enumerate(mu_array):
  for k,x in enumerate(carray):
    for i in range(0,500):
      rand = np.random.uniform(-3,3,3)
      #rand2 = np.random.uniform(1/10,10)

      Mv = 10**mu
      Mhc = Mv*(10**(rand[0]))
      M24t = Mv*(10**(rand[1]))
      M24o = Mv*(10**(rand[2]))

      lambda12=lambda_12(Mv,Mhc,M24t,10**mu,x)
      lambda23=lambda_23(Mv,Mhc,M24o,M24t,10**mu,x)

      ax1.plot(lambda12,lambda23,color=colors[k], marker='o')

      #if (Mv < 1e16):
        #ax1.plot(lambda12,lambda23,color=colors[j], marker = '*')#, label=r'$M_{V} < 10^{16}$ GeV, $\mu = 10^{%i}$ GeV' %mu)
      #else:
        #ax1.plot(lambda12,lambda23,color=colors[j], marker='o')#, label=r'$M_{V} > 10^{16}$ GeV, $\mu = 10^{%i}$ GeV'  %mu)


ax1.set_xlabel(r'$\Delta\lambda_{12}$',fontsize=14)
ax1.set_ylabel(r'$\Delta\lambda_{23}$',fontsize=14)
#ax1.set_title(r'Threshold corrections in minimal SU(5) with $\mu =10^{%i}$ GeV' %mu_array[0])
ax1.set_title(r'Threshold corrections in minimal SU(5)')
ax1.legend(loc=2)


plt.grid()
plt.show()
