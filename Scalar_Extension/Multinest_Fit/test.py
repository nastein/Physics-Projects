import json
import numpy
import numpy as np
from numpy import log, exp, pi, sqrt
import scipy.stats, scipy
import pymultinest
import matplotlib.pyplot as plt
import os
import random
import sys

#Ok so my problem has 3(N) + 1 parameters (quartic,vev,mixing for each scalar plus the higgs quartic)

def main():
  randomseed = sys.argv[2]
  random.seed(randomseed)

  outfolder = '/Users/noahsteinberg/Physics/James_Research/Scalar_Extension/Multinest_Fit/results/'
  if not os.path.exists(outfolder):
      os.makedirs(outfolder)

  N = int(sys.argv[1])
  vev = 246
  higgs_unc = 3*.14

  quadconst = 1
  quad_coupling = quadconst*pi
  xi_min = 1e3
  xi_max = 1e4
  epsilon = .01

  eps = epsilon*random.uniform(-1,1)

  C = 10
  mu = 5


  def rho(M):
    r = ((M**4)*(.000580808*C - .000580808*log(M*M/mu*mu) + .00195404) + (M**2)*(-8.58255*C + 8.58255*log(M*M/mu*mu) - 17.1651*log(M) + 51.8896) + 31207.1*C - 31207.1*log(M*M/mu*mu) + 62414.3*log(M) - 204116)/((M**4) - 14776.9*(M**2) + 5.37306*(10**7))
    return r

  format_list = [str(eps),str(xi_min),str(xi_max),str(quadconst)]
  print('Running scan with: epsilon = {}, Ximin = {}, Ximax = {}, quartic = {}pi.'.format(*format_list))
  print('Scan has N = {} additional scalars.'.format(str(N)))

  theta_min = []
  theta_max = []

  for i in range(0,N):
    theta_min.append(-quad_coupling)
    theta_min.append(xi_min)
    theta_max.append(quad_coupling)
    theta_max.append(xi_max)
  theta_min.append(-quad_coupling)
  theta_max.append(quad_coupling)

  theta_interval = list(np.array(theta_max)-np.array(theta_min))


  def prior(cube, ndim, nparams):
    for i in range(ndim):
      cube[i] = cube[i] * theta_interval[i] + theta_min[i]
    return cube

  def loglike(cube, ndim, nparams):

    likelihood1 = 0
    likelihood2 = 0

    lam = cube[2*N]
    higgs_mass_sq = 2*lam*vev*vev
    m_sq = []

    for i in range(0,N):
      kap = cube[2*i]
      xi = cube[1 + 2*i]
      m_sq.append(2*kap*(eps**2) + ((vev*xi*eps)**2)/(2*kap*xi*xi - 2*vev*vev*lam))
      higgs_mass_sq += ((vev*xi*eps)**2)/(2*vev*vev*lam - 2*kap*xi*xi)


    likelihood1 = exp(-(((125.7**2) - higgs_mass_sq) / (2*125.7*higgs_unc))**2)

    for i in range(0,N):
      likelihood2 += min(0,m_sq[i])

    return (log(likelihood1) + likelihood2)


  parameters = []
  for i in range(0,N):
    parameters.append('kappa_'+str(i))
    parameters.append('xi_'+str(i))
  parameters.append('lam')
  n_params = len(parameters)

  pymultinest_options = {'importance_nested_sampling': False, 'resume': False, 'verbose': True,'sampling_efficiency': 'model','init_MPI': False, 'evidence_tolerance': 0.5,'const_efficiency_mode': False}

  nlive=500

  pymultinest.run(loglike, prior, n_params, outputfiles_basename = outfolder, n_live_points = nlive, **pymultinest_options)

  json.dump(parameters, open(outfolder + 'params.json', 'w'))

  An = pymultinest.Analyzer(n_params, outputfiles_basename = outfolder)
  best_fit_params = An.get_best_fit()['parameters']
  max_LL_multinest = An.get_best_fit()['log_likelihood']
  s = An.get_stats()

  lam = best_fit_params[2*N]
  higgs_mass_sq = 2*lam*vev*vev

  m_sq = []
  higgs_sup = 0
  delta_rho = 0

  for i in range(0,N):
      kap = best_fit_params[2*i]
      xi = best_fit_params[1 + 2*i]
      m_sq.append(2*kap*(eps**2) + ((vev*xi*eps)**2)/(2*kap*xi*xi - 2*vev*vev*lam))
      higgs_mass_sq += ((vev*xi*eps)**2)/(2*vev*vev*lam - 2*kap*xi*xi)

      higgs_sup += (vev*vev*eps*eps*xi*xi)/((2*vev*vev*lam - 2*kap*xi*xi)**2)
      delta_rho += rho(sqrt(m_sq[i]))*(vev*vev*eps*eps*xi*xi)/((2*vev*vev*lam - 2*kap*xi*xi)**2)

  higgs_sup = 1/(sqrt(1 + higgs_sup))
  delta_rho = delta_rho*pow(higgs_sup,2) + rho(sqrt(higgs_mass_sq))*pow(higgs_sup,2) - rho(sqrt(higgs_mass_sq))

  print('Higgs suppression = {}'.format(str(higgs_sup)))
  print('Delta rho = {}'.format(str(delta_rho)))
  print('Higgs mass = ', sqrt(higgs_mass_sq))
  print('Scalar mass = ', sqrt(m_sq))

if __name__ == "__main__":
  main()
