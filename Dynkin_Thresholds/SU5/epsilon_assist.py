import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import sqrt as sqrt
from numpy import log as log
pi = np.pi

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

Mpl = 2.435e18
c1 = -1/3
c2 = -1
c3 = 2/3

cs = [.001,.01,.1,10,100]
#cs = np.linspace(0.001,.01,20)
#cs = [-1]
#print(cs)
#cs = [.01]

mu_array = [14, 14.5, 15, 15.5, 16, 16.5, 17]
sm12 = {14: -107.609, 14.5: -157.762, 15: -207.919, 15.5: -258.078, 16: -308.24, 16.5: -358.404, 17: -408.571}
sm23 = {14: 134.802, 14.5: 107.787, 15: 80.7879, 15.5: 53.8046, 16: 26.8361, 16.5: -0.118745, 17: -27.0608}



def unified(Mv,Mhc,M24t,M24o,mu,c):
  #onetwo = ((lambda_12(Mv,Mhc,M24t,mu,c) - sm12[mu])/sm12[mu])**2
  #twothree = ((lambda_23(Mv,Mhc,M24o,M24t,mu,c) - sm23[mu])/sm23[mu])**2
  onetwo = (lambda_12(Mv,Mhc,M24t,mu,c) - sm12[mu])**2
  twothree = (lambda_23(Mv,Mhc,M24o,M24t,mu,c) - sm23[mu])**2

  return np.sqrt(onetwo + twothree)

def lambda_12(Mv,Mhc,M24t,mu,c):
  muu = 10**mu
  return -(14/3) + 112*log(Mv/muu) + 2*log(M24t/muu) - (8/6)*log(Mhc/muu) + 8*c*(Mv/Mpl)*(c1-c2)*(48*np.pi*np.pi)

def lambda_23(Mv,Mhc,M24o,M24t,mu,c):
  muu = 10**mu
  return -2 + 21*log(Mv/muu) + (1/2)*log(Mhc*M24o/(muu*muu)) - 2*log(M24t/muu) + 8*c*(Mv/Mpl)*(c2-c3)*(48*np.pi*np.pi)


# rand = np.random.uniform(-3,3,3)
# mu = 16
# Mv = 10**mu
# Mhc = Mv*(10**rand[0])
# M24t = Mv*(10**rand[1])
# M24o = Mv*(10**rand[2])
# print('Mv = ', Mv/1e16)
# print('Mhc = ', Mhc/1e16)
# print('M24t =', M24t/1e16)
# print('M24o =', M24o/1e16)
# lambda12 = lambda_12(Mv,Mhc,M24t,mu,cs[0])
# lambda23 = lambda_23(Mv,Mhc,M24o,M24t,mu,cs[0])
# un = unified(Mv,Mhc,M24t,M24o,mu,cs[0])
# print('lambda12 = ', lambda12, ', lambda23 = ', lambda23, ', unified = ', un)

fig = plt.figure()

ax1 = fig.add_subplot(1,1,1)
ax1.scatter([],[],color='black')

for c in cs:
  print('c = ', c)
  for j,mu in enumerate(mu_array):
    for i in range(0,1000):
      rand = np.random.uniform(-2,2,3)
      Mv = 10**mu
      Mhc = Mv*(10**rand[0])
      M24t = Mv*(10**rand[1])
      M24o = Mv*(10**rand[2])

      #print(unified(Mv,Mhc,M24t,M24o,mu,c))
      if (unified(Mv,Mhc,M24t,M24o,mu,c) < 200):
        ax1.plot(np.log10(Mv), c, 'bo')


ax1.set_xlabel('log10(Mv/GeV)', fontsize=14)
ax1.set_ylabel('Dim 5 coefficient', fontsize=14)
ax1.set_yscale('log')
ax1.set_xlim([13.5,17.5])

plt.grid()
plt.show()

