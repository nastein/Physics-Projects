import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import sqrt as sqrt
from numpy import log as log
pi = np.pi

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

MGUT = 1e16
mu_array = [15,16,17]
colors = ['r','b','g','y']
markers = ['*', 'o']

# def lambda_12V(Mv,mu):
#   return -16/3 + (21*16/3)*log(Mv/mu)
# def lambda_23V(Mv,mu):
#   return -1 + 21*log(Mv/mu)
# def lambda_12trip(Mhc,mu):
#   return -4/3 -(4/3)*log(Mhc/mu)
# def lambda_23trip(Mhc,mu):
#   return .5 + .5*log(Mhc/mu)
# def lambda_12Adj(M24o,mu):
#   return 2 + 2*log(M24o/mu)
# def lambda_23Adj(M24o, M24t, mu):
#   return -(3/2) -(3/2)*log(((M24o/mu)**(1/2))*((mu/M24t)**2))

fig = plt.figure()

# ax2 = fig.add_subplot(111)
# ax2.scatter([],[],color='red')

# for i in range(0,3000):
#   rand = np.random.uniform(-10,10,5)
#   Mv = MGUT*10**(rand[0])
#   Mhc = Mv*10**(rand[1])
#   M24t = Mv*10**(rand[2])
#   M24o = Mv*10**(rand[3])

#   ax2.plot(lambda_12V(Mv,1e16),lambda_23V(Mv,1e16),color='red',marker='*', label='Vector')
#   ax2.plot(lambda_12trip(Mhc,1e16),lambda_23trip(Mhc,1e16),color='green', marker='o', label='Triplet Higgs')
#   ax2.plot(lambda_12Adj(M24o,1e16), lambda_23Adj(M24o,M24t,1e16),color='blue', marker='x', label='Adjoint Higgs')


#   #V = plt.arrow(0,0,lambda_12V(Mv,1e16),lambda_23V(Mv,1e16),color='red', label='Vector')
#   #trip = plt.arrow(0,0,lambda_12trip(Mhc,1e16),lambda_23trip(100*Mhc,1e16),color='green', label='Triplet Higgs')
#   #adj = plt.arrow(0,0,lambda_12Adj(M24o,1e16), lambda_23Adj(M24o,M24t,1e16),color='blue', label='Adjoint Higgs')


# ax2.set_xlabel(r'$\Delta\lambda_{12}$',fontsize=14)
# ax2.set_ylabel(r'$\Delta\lambda_{23}$',fontsize=14)
# ax2.set_title(r'$\Delta\lambda$ vectors in minimal SU(5), $\mu = 10^{16}$ GeV')
# #ax2.set_xlim(-10,10)
# #ax2.set_ylim(-5,5)
# #ax2.legend([V,trip,adj],['Vector','Triplet Higgs','Adjoint Higgs'],loc=2)
# ax2.legend(loc=2)

def lambda_12(Mv,Mhc,M24t,mu):
  return -(14/3) + 112*log(Mv/mu) + 2*log(M24t/mu) - (8/6)*log(Mhc/mu)

def lambda_23(Mv,Mhc,M24o,M24t,mu):
  return -2 + 21*log(Mv/mu) + (1/2)*log(Mhc*M24o/(mu*mu)) - 2*log(M24t/mu)

sm12 = [293.51, 193.251, 92.977, -7.30987, -107.609, -207.919, -308.24, -408.571, -508.912]
sm23 = [351.668, 297.305, 243.049, 188.885, 134.802, 80.7879, 26.8361, -27.0608, -80.9081]
text = [r'$10^{10}$',r'$10^{11}$',r'$10^{12}$',r'$10^{13}$',r'$10^{14}$',r'$10^{15}$',r'$10^{16}$',r'$10^{17}$',r'$10^{18}$']

ax1 = fig.add_subplot(1,1,1)
ax1.scatter(sm12,sm23,color='black',zorder=2, label='SM threshold correction')

for i, txt in enumerate(text):
  ax1.annotate(txt, (sm12[i], sm23[i]))

rand = np.random.uniform(-3,2,5)
Mv = MGUT*(10**rand[0])
Mhc = MGUT*10**(rand[1])
M24t = MGUT*10**(rand[2])
M24o = MGUT*10**(rand[3])

for j,mu in enumerate(mu_array):
  #for i in range(0,500):
    #rand = np.random.uniform(-3,2,5)
    #Mv = MGUT*(10**rand[0])
    #Mhc = MGUT*10**(rand[1])
    #M24t = MGUT*10**(rand[2])
    #M24o = MGUT*10**(rand[3])


    lambda12=lambda_12(Mv,Mhc,M24t,10**mu)
    lambda23=lambda_23(Mv,Mhc,M24o,M24t,10**mu)
    #print('lambda_23 = ', lambda23)
    #print('lambda_12 = ', lambda12)

    if (Mv < 1e16):
      ax1.plot(lambda12,lambda23,color=colors[j], marker = '*', label=r'$M_{V} < 10^{16}$ GeV, $\mu = 10^{%i}$ GeV' %mu)
    else:
      ax1.plot(lambda12,lambda23,color=colors[j], marker='o', label=r'$M_{V} > 10^{16}$ GeV, $\mu = 10^{%i}$ GeV'  %mu)

ax1.plot(sm12[6],sm23[6],marker='X', markersize=13, color='green')


dx15 = sm12[5] - lambda_12(Mv,Mhc,M24t,10**15)
dy15 = sm23[5] - lambda_23(Mv,Mhc,M24o, M24t,10**15)
ax1.arrow(lambda_12(Mv,Mhc,M24t,10**15),lambda_23(Mv,Mhc,M24o,M24t,10**15), dx15, dy15, shape='full')

dx16 = sm12[6] - lambda_12(Mv,Mhc,M24t,10**16)
dy16 = sm23[6] - lambda_23(Mv,Mhc,M24o, M24t,10**16)
ax1.arrow(lambda_12(Mv,Mhc,M24t,10**16),lambda_23(Mv,Mhc,M24o,M24t,10**16), dx16, dy16, shape='full')

dx17 = sm12[7] - lambda_12(Mv,Mhc,M24t,10**17)
dy17 = sm23[7] - lambda_23(Mv,Mhc,M24o, M24t,10**17)
ax1.arrow(lambda_12(Mv,Mhc,M24t,10**17),lambda_23(Mv,Mhc,M24o,M24t,10**17), dx17, dy17, shape='full')

ax1.set_xlabel(r'$\Delta\lambda_{12}$',fontsize=14)
ax1.set_ylabel(r'$\Delta\lambda_{23}$',fontsize=14)
#ax1.set_title(r'Threshold corrections in minimal SU(5) with $\mu =10^{%i}$ GeV' %mu_array[0])
ax1.set_title(r'Threshold corrections in minimal SU(5)')
ax1.legend(loc=2)

print(r'$Mv/(1e16 GeV)$  =', Mv/1e16)
print(r'$Mhc/(1e16 GeV)$  = ', Mhc/1e16)
print(r'$M24t/(1e16 GeV)$  = ', M24t/1e16)
print(r'$M24o/(1e16 GeV)$  = ', M24o/1e16)

plt.text(sm12[5] - dx15/2,sm23[5] - dy15/2, '(%f, %f)' %(dx15,dy15))
plt.text(sm12[6] - dx16/2,sm23[6] - dy16/2, '(%f, %f)' %(dx16,dy16))
plt.text(sm12[7] - dx17/2,sm23[7] - dy17/2, '(%f, %f)' %(dx17,dy17))

plt.grid()
plt.show()
