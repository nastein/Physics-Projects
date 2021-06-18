import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import operator

import matplotlib as mpl
import matplotlib.pyplot as plt
# ###########################
# # Setup Plotting Defaults #
# ###########################
# # For more options see https://matplotlib.org/users/customizing.html
# # Commands for high detail plots (much larger in file size though)
# #mpl.rcParams['agg.path.chunksize'] = 1000
# #mpl.rcParams['savefig.dpi'] = 1000
# # Line styles
# mpl.rcParams['lines.linewidth'] = 1.5
# mpl.rcParams['lines.antialiased'] = True
# mpl.rcParams['lines.dashed_pattern'] = 2.8, 1.5
# mpl.rcParams['lines.dashdot_pattern'] = 4.8, 1.5, 0.8, 1.5
# mpl.rcParams['lines.dotted_pattern'] = 1.1, 1.1
# mpl.rcParams['lines.scale_dashes'] = True
# # Default colors
# from cycler import cycler
# mpl.rcParams['axes.prop_cycle'] = cycler('color',['cornflowerblue','forestgreen','maroon','goldenrod','firebrick','mediumorchid', 'navy', 'brown'])
# # Fonts
# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['font.serif'] = 'CMU Serif'
# mpl.rcParams['font.sans-serif'] = 'CMU Sans Serif, DejaVu Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif'
# mpl.rcParams['text.usetex'] = True
# # Axes
# mpl.rcParams['axes.linewidth'] = 1.0
# mpl.rcParams['axes.labelsize'] = 24
# mpl.rcParams['axes.labelpad'] = 9.0
# # Tick marks - the essence of life
# mpl.rcParams['xtick.top'] = True
# mpl.rcParams['xtick.major.size'] = 10
# mpl.rcParams['xtick.minor.size'] = 5
# mpl.rcParams['xtick.major.width'] = 1.0
# mpl.rcParams['xtick.minor.width'] = 0.75
# mpl.rcParams['xtick.major.pad'] = 8
# mpl.rcParams['xtick.labelsize'] = 22
# mpl.rcParams['xtick.direction'] = 'in'
# mpl.rcParams['xtick.minor.visible'] = True
# mpl.rcParams['ytick.right'] = True
# mpl.rcParams['ytick.major.size'] = 10
# mpl.rcParams['ytick.minor.size'] = 5
# mpl.rcParams['ytick.major.width'] = 1.0
# mpl.rcParams['ytick.minor.width'] = 0.75
# mpl.rcParams['ytick.major.pad'] = 8
# mpl.rcParams['ytick.labelsize'] = 22
# mpl.rcParams['ytick.direction'] = 'in'
# mpl.rcParams['ytick.minor.visible'] = True
# # Legend
# mpl.rcParams['legend.fontsize'] = 14
# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['legend.framealpha'] = 0.8
# #mpl.rcParams['legend.edgecolor'] = 'black'
# mpl.rcParams['legend.fancybox'] = True
# mpl.rcParams['legend.borderpad'] = 0.4 # border whitespace
# mpl.rcParams['legend.labelspacing'] = 0.5 # the vertical space between the legend entries
# mpl.rcParams['legend.handlelength'] = 1.5 # the length of the legend lines
# mpl.rcParams['legend.handleheight'] = 0.7 # the height of the legend handle
# mpl.rcParams['legend.handletextpad'] = 0.5 # the space between the legend line and legend text
# mpl.rcParams['legend.borderaxespad'] = 0.5 # the border between the axes and legend edge
# mpl.rcParams['legend.columnspacing'] = 2.0 # column separation
# # Figure size
# mpl.rcParams['figure.figsize'] = 10, 6
# # Save details
# mpl.rcParams['savefig.bbox'] = 'tight'
# mpl.rcParams['savefig.pad_inches'] = 0.1

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

LEPI_2A = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/LEPI_2A.csv", delimiter=",", unpack=True)
LEPI_2AM = tuple([x for x in LEPI_2A[0]])
LEPI_2AG = tuple([x*1000 for x in LEPI_2A[2]])
t = sorted(zip(LEPI_2AM, LEPI_2AG), key=operator.itemgetter(0))
LEPI_2AM, LEPI_2AG = zip(*t)
LEPI_2AL = tuple([x**-1 for x in LEPI_2AG])

pp3aAtlas = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/pp3aAtlas.csv", delimiter=",", unpack=True)
pp3aAtlasM = tuple([x for x in pp3aAtlas[0]])
pp3aAtlasG = tuple([x*1000 for x in pp3aAtlas[2]])
t = sorted(zip(pp3aAtlasM, pp3aAtlasG), key=operator.itemgetter(0))
pp3aAtlasM, pp3aAtlasG = zip(*t)
pp3aAtlasL = tuple([x**-1 for x in pp3aAtlasG])

OPAL_3A = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/OPAL_3A.csv", delimiter=",", unpack=True)
OPAL_3AM = tuple([x for x in OPAL_3A[0]])
OPAL_3AG = tuple([x*1000 for x in OPAL_3A[2]])
t = sorted(zip(OPAL_3AM, OPAL_3AG), key=operator.itemgetter(0))
OPAL_3AM, OPAL_3AG = zip(*t)
OPAL_3AL = tuple([x**-1 for x in OPAL_3AG])

ppZ3aAtlas = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/ppZ3aAtlas.csv", delimiter=",", unpack=True)
ppZ3aAtlasM = tuple([x for x in ppZ3aAtlas[0]])
ppZ3aAtlasG = tuple([x*1000 for x in ppZ3aAtlas[2]])
t = sorted(zip(ppZ3aAtlasM, ppZ3aAtlasG), key=operator.itemgetter(0))
ppZ3aAtlasM, ppZ3aAtlasG = zip(*t)
ppZ3aAtlasL = tuple([x**-1 for x in ppZ3aAtlasG])

CMSPbPb = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/CMSPbPb.csv", delimiter=",", unpack=True)
CMSPbPbM = tuple([x for x in CMSPbPb[0]])
CMSPbPbG = tuple([x*1000 for x in CMSPbPb[2]])
t = sorted(zip(CMSPbPbM, CMSPbPbG), key=operator.itemgetter(0))
CMSPbPbM, CMSPbPbG = zip(*t)
CMSPbPbL = tuple([x**-1 for x in CMSPbPbG])

AtlasPbPb = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/AtlasPbPb.csv", delimiter=",", unpack=True)
AtlasPbPbM = tuple([x for x in AtlasPbPb[0]])
AtlasPbPbG = tuple([x*1000 for x in AtlasPbPb[2]])
t = sorted(zip(AtlasPbPbM, AtlasPbPbG), key=operator.itemgetter(0))
AtlasPbPbM, AtlasPbPbG = zip(*t)
AtlasPbPbL = tuple([x**-1 for x in AtlasPbPbG])

beamdump = np.loadtxt("/Users/noahsteinberg/Desktop/collider_data/Beamdump.csv", delimiter=",", unpack=True)
beamdumpM = tuple([x for x in beamdump[0]])
beamdumpG = tuple([x*1000 for x in beamdump[2]])
t = sorted(zip(beamdumpM, beamdumpG), key=operator.itemgetter(1))
beamdumpM, beamdumpG = zip(*t)
beamdumpL = tuple([x**-1 for x in beamdumpG])

#ILCM, ILCG, ILCL = np.loadtxt("/Users/noahsteinberg/Physics/James_Research/ZDecay/ILCdata.csv", delimiter=",", unpack=True)
#ILCG = tuple([x for x in ILCG])
#ILCL = tuple([x for x in ILCL])
#ILCM, ILCG, ILCL = np.loadtxt("/Users/noahsteinberg/Desktop/ILCexclusion.csv", delimiter=",", unpack=True)
ILCM = [10,15,20,25,30,40,50,75,100,125,150,200]
ILCG = [.01344, .109, .1, .0932, .08794, .08222, .07476, .066967, .061609, .0655871, .0676537, .0668096]
ILCL = tuple([x**-1 for x in ILCG])

fig, (ax) = plt.subplots(1,sharex=True)
ax.set_yscale('log')
ax.set_ylim(1, 1e2)
ax.set_xscale('linear')
ax.set_xlim((5, 200))
ax.spines['left'].set_visible(False)
ax.yaxis.set_ticks_position('right')
ax.yaxis.set_visible(False)
divider = make_axes_locatable(ax)
axLin = divider.append_axes("left", size=2.0, pad=0, sharey=ax)
axLin.set_xscale('log')
axLin.set_xlim((.035, 5))
axLin.tick_params(axis="x", labelsize=22)
axLin.tick_params(axis="y", labelsize=22)
axLin.spines['right'].set_visible(False)
axLin.yaxis.set_ticks_position('left')
axLin.set_ylabel(r'$g^{-1}_{aBB} [\rm TeV]$', fontsize=23)
plt.setp(axLin.get_xticklabels(), visible=True)
ax.tick_params(axis="x", labelsize=22)
ax.tick_params(axis="y", labelsize=22)
ax.set_ylabel(r'$g^{-1}_{aBB} [\rm TeV]$', fontsize=23)
ax.set_xlabel(r'$m_{a} [\rm GeV]$', fontsize=23)

ax.plot(beamdumpM, beamdumpL, color='pink', label=r'Beam Dump')
axLin.plot(beamdumpM, beamdumpL, color='pink', label=r'Beam Dump')
ax.fill_between(beamdumpM, y2 = 100000, y1 = beamdumpL,where=None, interpolate=True, color='pink', alpha=0.2)
axLin.fill_between(beamdumpM, y2 = 100000, y1 = beamdumpL,where=None, interpolate=True, color='pink', alpha=0.2)

ax.plot(LEPI_2AM, LEPI_2AL, color='cyan', label=r'LEPI $e^{+}e^{-}\rightarrow 2\gamma$')
axLin.plot(LEPI_2AM, LEPI_2AL, color='cyan', label=r'LEPI $e^{+}e^{-}\rightarrow 2\gamma$')
ax.fill_between(LEPI_2AM, y1 = 0, y2 = LEPI_2AL,where=None, interpolate=True, color='cyan', alpha=0.2)
axLin.fill_between(LEPI_2AM, y1 = 0, y2 = LEPI_2AL,where=None, interpolate=True, color='cyan', alpha=0.2)

ax.plot(OPAL_3AM, OPAL_3AL, color='blue', label=r'OPAL $e^{+}e^{-}\rightarrow 3\gamma$')
axLin.plot(OPAL_3AM, OPAL_3AL, color='blue', label=r'OPAL $e^{+}e^{-}\rightarrow 3\gamma$')
ax.fill_between(OPAL_3AM, y1 = 0, y2 = OPAL_3AL,where=None, interpolate=True, color='blue', alpha=0.2)
axLin.fill_between(OPAL_3AM, y1 = 0, y2 = OPAL_3AL,where=None, interpolate=True, color='blue', alpha=0.2)

ax.plot(pp3aAtlasM, pp3aAtlasL, color='orange', label=r'Atlas pp $\rightarrow 3\gamma$')
axLin.plot(pp3aAtlasM, pp3aAtlasL, color='orange', label=r'Atlas pp $\rightarrow 3\gamma$')
ax.fill_between(pp3aAtlasM, y1 = 0, y2 = pp3aAtlasL, where=None, interpolate=True, color='orange', alpha=0.2)
axLin.fill_between(pp3aAtlasM, y1 = 0, y2 = pp3aAtlasL, where=None, interpolate=True, color='orange', alpha=0.2)

ax.plot(ppZ3aAtlasM, ppZ3aAtlasL, color='grey', label=r'Atlas pp $\rightarrow Z(3\gamma)$')
axLin.plot(ppZ3aAtlasM, ppZ3aAtlasL, color='grey', label=r'Atlas pp $\rightarrow Z(3\gamma)$')
ax.fill_between(ppZ3aAtlasM, y1 = 0, y2 = ppZ3aAtlasL, where=None, interpolate=True, color='blue', alpha=0.2)
axLin.fill_between(ppZ3aAtlasM, y1 = 0, y2 = ppZ3aAtlasL, where=None, interpolate=True, color='blue', alpha=0.2)

ax.plot(CMSPbPbM, CMSPbPbL, color='green', label=r'CMS PbPb')
axLin.plot(CMSPbPbM, CMSPbPbL, color='green', label=r'CMS PbPb')
ax.fill_between(CMSPbPbM, y1 = 0, y2 = CMSPbPbL,where=None, interpolate=True, color='green', alpha=0.2)
axLin.fill_between(CMSPbPbM, y1 = 0, y2 = CMSPbPbL,where=None, interpolate=True, color='green', alpha=0.2)

ax.plot(AtlasPbPbM, AtlasPbPbL, color='red', label=r'Atlas PbPb')
axLin.plot(AtlasPbPbM, AtlasPbPbL, color='red', label=r'Atlas PbPb')
ax.fill_between(AtlasPbPbM, y1 = 0, y2 = AtlasPbPbL,where=None, interpolate=True, color='red', alpha=0.2)
axLin.fill_between(AtlasPbPbM, y1 = 0, y2 = AtlasPbPbL,where=None, interpolate=True, color='red', alpha=0.2)

ax.plot(ILCM, ILCL, color ='purple', label=r'ILC 500 GeV $(5\,\rm{ab^{-1})}$')
axLin.plot(ILCM, ILCL, color ='purple', label=r'ILC 500 GeV $(5\,\rm{ab^{-1})}$')
ax.fill_between(ILCM, y1 = 0, y2 = ILCL,where=None, interpolate=True, color='purple', alpha=0.2)
axLin.fill_between(ILCM, y1 = 0, y2 = ILCL,where=None, interpolate=True, color='purple', alpha=0.2)

ax.axvline(x=5, linewidth = 3, color='black', linestyle='--')

ax.legend(prop=dict(size=19))

fig2, (ax2) = plt.subplots(1,sharex=True)
ax2.set_yscale('log')
ax2.set_ylim(.01, 1)
ax2.set_xscale('linear')
ax2.set_xlim((5, 200))
ax2.spines['left'].set_visible(False)
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_visible(False)
divider = make_axes_locatable(ax2)
axLin2 = divider.append_axes("left", size=2.0, pad=0, sharey=ax2)
axLin2.set_xscale('log')
axLin2.set_xlim((.035, 5))
axLin2.tick_params(axis="x", labelsize=22)
axLin2.tick_params(axis="y", labelsize=22)
axLin2.spines['right'].set_visible(False)
axLin2.yaxis.set_ticks_position('left')
axLin2.set_ylabel(r'$g_{aBB} [\rm TeV^{-1}]$', fontsize=23)
ax2.set_xlabel(r'$m_{a} [GeV]$', fontsize=23)
plt.setp(axLin2.get_xticklabels(), visible=True)
ax2.tick_params(axis="x", labelsize=22)
ax2.tick_params(axis="y", labelsize=22)
ax2.set_ylabel(r'$g_{aBB} [\rm TeV^{-1}]$', fontsize=23)
ax2.set_xlabel(r'$m_{a} [\rm GeV]$', fontsize=23)


ax2.plot(beamdumpM, beamdumpG, color='pink', label=r'Beam Dump')
axLin2.plot(beamdumpM, beamdumpG, color='pink', label=r'Beam Dump')
ax2.fill_between(beamdumpM, y1 = .00000001, y2 = beamdumpG,where=None, interpolate=True, color='pink', alpha=0.2)
axLin2.fill_between(beamdumpM, y1 = .00000001, y2 = beamdumpG,where=None, interpolate=True, color='pink', alpha=0.2)

ax2.plot(LEPI_2AM, LEPI_2AG, color='cyan', label=r'LEPI $e^{+}e^{-}\rightarrow 2\gamma$')
axLin2.plot(LEPI_2AM, LEPI_2AG, color='cyan', label=r'LEPI $e^{+}e^{-}\rightarrow 2\gamma$')
ax2.fill_between(LEPI_2AM, y2 = 100, y1 = LEPI_2AG,where=None, interpolate=True, color='cyan', alpha=0.2)
axLin2.fill_between(LEPI_2AM, y2 = 100, y1 = LEPI_2AG,where=None, interpolate=True, color='cyan', alpha=0.2)

ax2.plot(OPAL_3AM, OPAL_3AG, color='blue', label=r'OPAL $e^{+}e^{-}\rightarrow 3\gamma$')
axLin2.plot(OPAL_3AM, OPAL_3AG, color='blue', label=r'OPAL $e^{+}e^{-}\rightarrow 3\gamma$')
ax2.fill_between(OPAL_3AM, y2 = 100, y1 = OPAL_3AG,where=None, interpolate=True, color='blue', alpha=0.2)
axLin2.fill_between(OPAL_3AM, y2 = 100, y1 = OPAL_3AG,where=None, interpolate=True, color='blue', alpha=0.2)

ax2.plot(pp3aAtlasM, pp3aAtlasG, color='orange', label=r'Atlas pp $\rightarrow 3\gamma$')
axLin2.plot(pp3aAtlasM, pp3aAtlasG, color='orange', label=r'Atlas pp $\rightarrow 3\gamma$')
ax2.fill_between(pp3aAtlasM, y2 = 100, y1 = pp3aAtlasG, where=None, interpolate=True, color='orange', alpha=0.2)
axLin2.fill_between(pp3aAtlasM, y2 = 100, y1 = pp3aAtlasG, where=None, interpolate=True, color='orange', alpha=0.2)

ax2.plot(ppZ3aAtlasM, ppZ3aAtlasG, color='grey', label=r'Atlas pp $\rightarrow Z(3\gamma)$')
axLin2.plot(ppZ3aAtlasM, ppZ3aAtlasG, color='grey', label=r'Atlas pp $\rightarrow Z(3\gamma)$')
ax2.fill_between(ppZ3aAtlasM, y2 = 100, y1 = ppZ3aAtlasG, where=None, interpolate=True, color='blue', alpha=0.2)
axLin2.fill_between(ppZ3aAtlasM, y2 = 100, y1 = ppZ3aAtlasG, where=None, interpolate=True, color='blue', alpha=0.2)

ax2.plot(CMSPbPbM, CMSPbPbG, color='green', label=r'CMS PbPb')
axLin2.plot(CMSPbPbM, CMSPbPbG, color='green', label=r'CMS PbPb')
ax2.fill_between(CMSPbPbM, y2 = 100, y1 = CMSPbPbG,where=None, interpolate=True, color='green', alpha=0.2)
axLin2.fill_between(CMSPbPbM, y2 = 100, y1 = CMSPbPbG,where=None, interpolate=True, color='green', alpha=0.2)

ax2.plot(AtlasPbPbM, AtlasPbPbG, color='red', label=r'Atlas PbPb')
axLin2.plot(AtlasPbPbM, AtlasPbPbG, color='red', label=r'Atlas PbPb')
ax2.fill_between(AtlasPbPbM, y2 = 100, y1 = AtlasPbPbG,where=None, interpolate=True, color='red', alpha=0.2)
axLin2.fill_between(AtlasPbPbM, y2 = 100, y1 = AtlasPbPbG,where=None, interpolate=True, color='red', alpha=0.2)

ax2.axvline(x=5, linewidth=3, color='black', linestyle='--')

ax2.plot(ILCM, ILCG, color ='purple', label=r'ILC 500 GeV $(5\,\rm{ab^{-1})}$')
axLin2.plot(ILCM, ILCG, color ='purple', label=r'ILC 500 GeV $(5\,\rm{ab^{-1})}$')
ax2.fill_between(ILCM, y2 = 100, y1 = ILCG,where=None, interpolate=True, color='purple', alpha=0.2)
axLin2.fill_between(ILCM, y2 = 100, y1 = ILCG,where=None, interpolate=True, color='purple', alpha=0.2)

ax2.legend(prop=dict(size=19))

plt.show()











