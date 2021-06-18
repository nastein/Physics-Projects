import json
import numpy
import numpy as np
from numpy import log, exp, pi, sqrt
import scipy.stats, scipy
import pymultinest
import matplotlib.pyplot as plt

#Ok so my problem has 3(N) + 1 parameters (quartic,vev,mixing for each scalar plus the higgs quartic)

if not os.path.exists(outfolder):
    os.makedirs(outfolder)

outfolder = '/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/Multinest_Fit/results/'

N = 3
vev = 246
higgs_unc = 3*.14

quad = 4*pi
xi_min = 1e4
xi_max = 1e5
epsilon = .1

theta_min = []
theta_max = []

for i in range(0,N):
  theta_min.append(-quad)
  theta_min.append(xi_min)
  theta_min.append(-epsilon)
  theta_max.append(quad)
  theta_max.append(xi_max)
  theta_max.append(epsilon)
theta_min.append(-quad)
theta_max.append(quad)

theta_interval = list(np.array(theta_max)-np.array(theta_min))

print(theta_interval)

def prior(cube, ndim, nparams):
  for i in range(ndim):
    cube[i] = cube[i] * theta_interval[i] + theta_min[i]
  return cube

def loglike(cube, ndim, nparams):

  lam = cube[3*N]
  higgs_mass_sq = 2*lam*vev*vev
  m_sq = []

  for i in range(0,N):
    kap = cube[3*i]
    xi = cube[1 + 3*i]
    eps = cube[2 + 3*i]
    m_sq.append(2*kap*(eps**2) + ((vev*xi*eps)**2)/(2*kap*xi*xi - 2*vev*vev*lam))
    higgs_mass_sq += ((vev*xi*eps)**2)/(2*vev*vev*lam - 2*kap*xi*xi)

  else:
    likelihood1 = exp(-(((125.7**2) - higgs_mass_sq) / (2*125.7*higgs_unc))**2)

  for i in range(0,N):
    likelihood1 *= (1 + np.tanh(m_sq[i]))

  return log(likelihood1)

# number of dimensions our problem has
parameters = ["kappa1","xi1","eps1","kappa2","xi2","eps2","kappa3","xi3","eps3","lam"]
n_params = len(parameters)

pymultinest_options = {'importance_nested_sampling': False, 'resume': False, 'verbose': True,'sampling_efficiency': 'model','init_MPI': False, 'evidence_tolerance': 0.5,'const_efficiency_mode': False}

nlive=500

# run MultiNest
pymultinest.run(loglike, prior, n_params, outputfiles_basename=outfolder,n_live_points = nlive,**pymultinest_options)

#json.dump(parameters, open('out/params.json', 'w'))

An = pymultinest.Analyzer(n_params,outputfiles_basename='out/')
best_fit_params = An.get_best_fit()['parameters']
max_LL_multinest = An.get_best_fit()['log_likelihood']
s = An.get_stats()

print('Best fit parameters = ', best_fit_params)

chain_file = 'out/' + 'post_equal_weights.dat'
samples = np.array(np.loadtxt(chain_file)[:,:-1])

medians = [s['marginals'][i]['median'] for i in range(n_params)]

print('Medians = ', medians)

lam = best_fit_params[3*N]
higgs_mass_sq = 2*lam*vev*vev

m_sq = []

for i in range(0,N):
    kap = best_fit_params[3*i]
    xi = best_fit_params[1 + 3*i]
    eps = best_fit_params[2 + 3*i]
    m_sq.append(2*kap*(eps**2) + ((vev*xi*eps)**2)/(2*kap*xi*xi - 2*vev*vev*lam))
    higgs_mass_sq += ((vev*xi*eps)**2)/(2*vev*vev*lam - 2*kap*xi*xi)

print('Higgs mass = ', sqrt(higgs_mass_sq))
print('Scalar mass = ', sqrt(m_sq))
