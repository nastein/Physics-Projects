import numpy as np
import matplotlib.pyplot as plt
from labellines import *

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})

mz = 91.2
mw = 80.379
cw = mw/mz

def G(m,g):
  return (g**2)*(m**3)*(cw**4)/(64*np.pi)

def En(m):
  return mz - ((mz**2 - m**2)/(2*mz))

def r(m):
  return En(m)/m

def l(m,g):
  return (1.98e-14)*r(m)/G(m,g)

gs = np.arange(-5.5, -2.9, .1)
gs = np.array(tuple([10**g for g in gs]))
ms = [0.4, 1.0, 5, 10.0]

m1 = np.array(tuple([l(ms[0],g) for g in gs]))
m2 = np.array(tuple([l(ms[1],g) for g in gs]))
m3 = np.array(tuple([l(ms[2],g) for g in gs]))
m4 = np.array(tuple([l(ms[3],g) for g in gs]))

fig, ax = plt.subplots()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-7,2.2e3)
ax.set_xlim(10**-5.5, 1e-3)
ax.tick_params(axis="x", labelsize=22)
ax.tick_params(axis="y", labelsize=22)
ax.set_ylabel(r'$c\tau$', fontsize=23)
ax.set_xlabel(r'$g_{aBB} [\rm GeV^{-1}]$', fontsize=23)

ax.plot(gs, m1, color='red', label='0.4 GeV')
ax.plot(gs, m2, color='blue', label='1.0 GeV')
ax.plot(gs, m3, color='green', label='5.0 GeV')
ax.plot(gs, m4, color='purple', label='10 GeV')
#ax.legend(prop=dict(size=16))

labelLines(plt.gca().get_lines(), align=False, fontsize=19)

plt.show()









