from __future__ import absolute_import, unicode_literals, print_function
import numpy as np
import pycuba
import matplotlib.pyplot as plt
#pycuba.demo()

if __name__ == '__main__':
  import math

  NDIM = 3
  NCOMP = 1

  NNEW = 1000
  NMIN = 2
  FLATNESS = 50.

  KEY1 = 47
  KEY2 = 1
  KEY3 = 1
  MAXPASS = 5
  BORDER = 0.
  MAXCHISQ = 10.
  MINDEVIATION = .25
  NGIVEN = 0
  LDXGIVEN = NDIM
  NEXTRA = 0
  MINEVAL = 0
  MAXEVAL = 50000

  KEY = 0

  me = 5.11e-4 #GeV
  eta = 2.4
  Q2max = 2.0;
  ptmin = 2.0;

  def Q2min(x):
    return (me**2)*((x**2)/(1 - x))

  def fmgw(w,see):
    return (0.0024887403141813187*((1.0444840000000003e-6*w*(0.5 - (957410.5491323943*see*(1 - (2*w)/math.sqrt(see)))/w**2))/math.sqrt(see) + (math.sqrt(see)*(2 - (4*w)/math.sqrt(see) + (4*w**2)/see)*math.log((1.9148210982647885e6*see*(1 - (2*w)/math.sqrt(see)))/w**2))/(2.*w)))/math.sqrt(see)

  def xmin(s,pt):
    return math.exp(-2*eta)*((1 + math.sqrt(1 - (4*(pt**2)/s)))/(1 - math.sqrt(1 - (4*(pt**2)/s))))

  def xmax(s,pt):
    return math.exp(2*eta)*((1 - math.sqrt(1 - (4*(pt**2)/s)))/(1 + math.sqrt(1 - (4*(pt**2)/s))))

  def ptmax(s):
    return math.sqrt(s)/2

  def g(x,s,see):
    return (1.0/(8.0*x))*fmgw(math.sqrt(s*x/4.0),see)*fmgw(math.sqrt(s/(4.0*x)),see)

  def dsigdpt(s,pt,g,m):
    return (1.8267964937623458*g**4*pt*(866.0699102689115*g**4*m**12 + g**8*m**16 + 375038.54473660036*m**6*s + 866.0699102689115*g**4*m**10*s + 375038.54473660036*m**2*pt**2*s**2 + m**4*s*(375038.54473660036*pt**2 + 187519.2723683001*s) + s**2*(187519.27236830018*pt**4 - 1.4551915228366852e-11*pt**2*s + 3.637978807091713e-12*s**2) + m**8*(187519.27236830018 + g**4*s*(-866.0699102689115*pt**2 + 433.03495513445574*s)))*(g**4*m**10*pt**2*s*(-144.3449850448186*pt**2 + 144.3449850448186*s) + m**6*pt**2*s*(-62506.42412276672*pt**2 + 62506.42412276672*s) + pt**4*s**3*(-312532.1206138336*pt**2 + 312532.1206138336*s) + g**8*m**16*(pt**4 - 2.*pt**2*s + s**2) + g**4*m**12*(866.0699102689117*pt**4 - 1732.1398205378234*pt**2*s + 866.0699102689117*s**2) + m**2*pt**2*s**2*(250025.6964910669*pt**4 - 500051.3929821338*pt**2*s + 250025.6964910669*s**2) + m**4*s*(-62506.42412276672*pt**6 - 375038.5447366003*pt**4*s - 187519.27236830015*pt**2*s**2 + 62506.42412276672*s**3) + m**8*(-144.3449850448186*g**4*pt**6*s + 187519.27236830015*s**2 + 144.3449850448186*g**4*s**4 + pt**4*(187519.27236830015 + 2598.2097308067355*g**4*s**2) + pt**2*(-375038.5447366003*s - 433.03495513445586*g**4*s**3))))/(math.sqrt(1 - (4*pt**2)/s)*(0.022791703727531375*g**4*m**8 + m**4*math.pi**2 - 19.739208802178716*m**2*s + math.pi**2*s**2)*(0.004618564797797922*g**4*m**12 + 5.33278519786454e-6*g**8*m**16 + 2.*m**6*s + 0.004618564797797922*g**4*m**10*s + 2.*m**2*pt**2*s**2 + pt**4*s**2 + m**4*s*(2.*pt**2 + s) + m**8*(1. + g**4*(-0.004618564797797922*pt**2 + 0.002309282398898961*s)*s))**2)

  def compute_xsec(sqrtsmin, sqrtsmax, sqrtsee, mass, gabb):

    see = sqrtsee**2
    smin = sqrtsmin**2
    smax = sqrtsmax**2

    def transform_s(s):
      return (smin + (smax - smin)*s)

    def transform_pt(pt,s):
      return (ptmin + (ptmax(transform_s(s)) - ptmin)*pt)

    def transform_x(x,s,pt):
      return (xmin(transform_s(s),transform_pt(pt,s)) + (xmax(transform_s(s),transform_pt(pt,s)) - xmin(transform_s(s),transform_pt(pt,s)))*x)

    def Integrand(ndim, xx, ncomp, ff, userdata):
      s,pt,x = [xx[i] for i in range(ndim.contents.value)]
      result = dsigdpt(transform_s(s),transform_pt(pt,s), gabb, mass)*g(transform_x(x,s,pt),transform_s(s),see)*(smax - smin)*(xmax(transform_s(s),transform_pt(pt,s)) - xmin(transform_s(s),transform_pt(pt,s)))*(ptmax(transform_s(s)) - ptmin)
      ff[0] = result
      return 0

    a = 0
    result = pycuba.Vegas(Integrand, NDIM, verbose=0)['results']
    for comp in result:
      a = float("%(integral).8f" % comp)
    return a


  m = list(range(5,100,2))
  #m = [50]
  xsecs500 = []
  xsecs250 = []
  for mass in m:
    xsecs500.append(compute_xsec(mass - 2, mass + 2, 500, mass, 1e-3))
    xsecs250.append(compute_xsec(mass - 2, mass + 2, 250, mass, 1e-3))


  plt.plot(m,xsecs500,color='red',label=r'$\sqrt{s} = 500$ GeV')
  plt.plot(m,xsecs250,color='black',label=r'$\sqrt{s} = 250$ GeV')
  plt.legend(loc="upper left")
  plt.show()



