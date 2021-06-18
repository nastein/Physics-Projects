import numpy as np
import csv
import os
import matplotlib.pyplot as plt


NH = 'NH'
N_directories = ['osc000','osc010','osc001']
IH = 'IH'
I_directories = ['osc100','osc110','osc101']

#Order is es, ibd, oss, es+oss
SM_NH = [[],[],[],[]]
FV_NH = [[],[],[],[]]
FP_NH = [[],[],[],[]]

SM_IH = [[],[],[],[]]
FV_IH = [[],[],[],[]]
FP_IH = [[],[],[],[]]

SM_NH_err = [[],[],[],[]]
FV_NH_err = [[],[],[],[]]
FP_NH_err = [[],[],[],[]]

SM_IH_err = [[],[],[],[]]
FV_IH_err = [[],[],[],[]]
FP_IH_err = [[],[],[],[]]

def ratio_err(a,b,siga,sigb):
  return np.sqrt(float(siga*siga)/float(b*b) + float(sigb*sigb)*(float(a*a)/float(b**4)))

for num,d in enumerate(N_directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/nu_events/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith("1.txt"):
        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          times = np.array(list(zip(*reader))[0], dtype=float)

          if (num == 0):
            if (filepath.endswith('es1.txt')):
              SM_NH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              SM_NH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              SM_NH[2].extend(times)
          if (num == 1):
            if (filepath.endswith('es1.txt')):
              FV_NH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              FV_NH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              FV_NH[2].extend(times)
          if (num == 2):
            if (filepath.endswith('es1.txt')):
              FP_NH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              FP_NH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              FP_NH[2].extend(times)

for num,d in enumerate(I_directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/nu_events/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith("1.txt"):# and ('es' in filepath):
        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          times = np.array(list(zip(*reader))[0], dtype=float)

          if (num == 0):
            if (filepath.endswith('es1.txt')):
              SM_IH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              SM_IH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              SM_IH[2].extend(times)
          if (num == 1):
            if (filepath.endswith('es1.txt')):
              FV_IH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              FV_IH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              FV_IH[2].extend(times)
          if (num == 2):
            if (filepath.endswith('es1.txt')):
              FP_IH[0].extend(times)
            if (filepath.endswith('ibd1.txt')):
              FP_IH[1].extend(times)
            if (filepath.endswith('o16e1.txt') or filepath.endswith('o16eb1.txt')):
              FP_IH[2].extend(times)

bins=range(1000,22000,1000)

for i in range(0,3):
  SM_NH[i],bin_edges = np.histogram(SM_NH[i],bins=bins)
  FV_NH[i],bin_edges = np.histogram(FV_NH[i],bins=bins)
  FP_NH[i],bin_edges = np.histogram(FP_NH[i],bins=bins)
  SM_IH[i],bin_edges = np.histogram(SM_IH[i],bins=bins)
  FV_IH[i],bin_edges = np.histogram(FV_IH[i],bins=bins)
  FP_IH[i],bin_edges = np.histogram(FP_IH[i],bins=bins)

  SM_NH_err[i] = np.sqrt(SM_NH[i])
  SM_IH_err[i] = np.sqrt(SM_IH[i])
  FV_NH_err[i] = np.sqrt(FV_NH[i])
  FV_IH_err[i] = np.sqrt(FV_IH[i])
  FP_NH_err[i] = np.sqrt(FP_NH[i])
  FP_IH_err[i] = np.sqrt(FP_IH[i])

SM_NH[3] = SM_NH[0] + SM_NH[2]
SM_IH[3] = SM_IH[0] + SM_IH[2]
FP_NH[3] = FP_NH[0] + FP_NH[2]
FP_IH[3] = FP_IH[0] + FP_IH[2]
FV_NH[3] = FV_NH[0] + FV_NH[2]
FV_IH[3] = FV_IH[0] + FV_IH[2]

SM_NH_err[3] = np.sqrt(SM_NH_err[0]**2 + SM_NH_err[2]**2)
SM_IH_err[3] = np.sqrt(SM_IH_err[0]**2 + SM_IH_err[2]**2)
FP_NH_err[3] = np.sqrt(FP_NH_err[0]**2 + FP_NH_err[2]**2)
FP_IH_err[3] = np.sqrt(FP_IH_err[0]**2 + FP_IH_err[2]**2)
FV_NH_err[3] = np.sqrt(FV_NH_err[0]**2 + FV_NH_err[2]**2)
FV_IH_err[3] = np.sqrt(FV_IH_err[0]**2 + FV_IH_err[2]**2)

fig1 = plt.figure(1)

ax1 = fig1.add_subplot(221)
ax1.errorbar(bin_edges[:-1],SM_NH[1], yerr=SM_NH_err[1], fmt='o-',label='Standard', markersize=5)
ax1.errorbar(bin_edges[:-1],FP_NH[1], yerr=FP_NH_err[1], fmt='o-',label='FP-NSSI', markersize=5)
ax1.errorbar(bin_edges[:-1],FV_NH[1], yerr=FV_NH_err[1], fmt='o-',label='FV-NSSI', markersize=5)

ax2 = fig1.add_subplot(222)
ax2.errorbar(bin_edges[:-1],SM_NH[3], yerr=SM_NH_err[3], fmt='o-',label='Standard', markersize=5)
ax2.errorbar(bin_edges[:-1],FP_NH[3], yerr=FP_NH_err[3], fmt='o-',label='FP-NSSI', markersize=5)
ax2.errorbar(bin_edges[:-1],FV_NH[3], yerr=FV_NH_err[3], fmt='o-',label='FV-NSSI', markersize=5)

ax3 = fig1.add_subplot(223)
ax3.errorbar(bin_edges[:-1],SM_IH[1], yerr=SM_NH_err[1], fmt='o-',label='Standard', markersize=5)
ax3.errorbar(bin_edges[:-1],FP_IH[1], yerr=FP_NH_err[1], fmt='o-',label='FP-NSSI', markersize=5)
ax3.errorbar(bin_edges[:-1],FV_IH[1], yerr=FV_NH_err[1], fmt='o-',label='FV-NSSI', markersize=5)

ax4 = fig1.add_subplot(224)
ax4.errorbar(bin_edges[:-1],SM_IH[3], yerr=SM_IH_err[3], fmt='o-',label='Standard', markersize=5)
ax4.errorbar(bin_edges[:-1],FP_IH[3], yerr=FP_IH_err[3], fmt='o-',label='FP-NSSI', markersize=5)
ax4.errorbar(bin_edges[:-1],FV_IH[3], yerr=FV_IH_err[3], fmt='o-',label='FV-NSSI', markersize=5)

xtick = [1000,3500,6000,8500,11000,13500,16000,18000,20000]

ax1.set_ylabel('Events/1000ms')
ax1.set_xlabel('Time (ms)')
ax1.set_title('IBD, NH', y=.9)
ax1.legend(loc='upper right');
ax1.set_xticks(xtick)
for tick in ax1.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax1.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax1.set_xlim(0,21000)

ax2.set_ylabel('Events/1000ms')
ax2.set_xlabel('Time (ms)')
ax2.set_title('ES + OSS, NH', y=.9)
ax2.set_xticks(xtick)
for tick in ax2.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax2.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax2.set_xlim(0,21000)

ax3.set_ylabel('Events/1000ms')
ax3.set_xlabel('Time (ms)')
ax3.set_title('IBD, IH', y=.9)
ax3.set_xticks(xtick)
for tick in ax3.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax3.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax3.set_xlim(0,21000)

ax4.set_ylabel('Events/1000ms')
ax4.set_xlabel('Time (ms)')
ax4.set_title('ES + OSS, IH', y=.9)
ax4.set_xticks(xtick)
for tick in ax4.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax4.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax4.set_xlim(0,21000)

fig2 = plt.figure(2)

ax1 = fig2.add_subplot(231)
ax1.errorbar(bin_edges[:-1],SM_NH[0], yerr=SM_NH_err[0], fmt='o-',label='Standard', markersize=5)
ax1.errorbar(bin_edges[:-1],FP_NH[0], yerr=FP_NH_err[0], fmt='o-',label='FP-NSSI', markersize=5)
ax1.errorbar(bin_edges[:-1],FV_NH[0], yerr=FV_NH_err[0], fmt='o-',label='FV-NSSI', markersize=5)

ax2 = fig2.add_subplot(232)
ax2.errorbar(bin_edges[:-1],SM_NH[1], yerr=SM_NH_err[1], fmt='o-',label='Standard', markersize=5)
ax2.errorbar(bin_edges[:-1],FP_NH[1], yerr=FP_NH_err[1], fmt='o-',label='FP-NSSI', markersize=5)
ax2.errorbar(bin_edges[:-1],FV_NH[1], yerr=FV_NH_err[1], fmt='o-',label='FV-NSSI', markersize=5)

ax3 = fig2.add_subplot(233)
ax3.errorbar(bin_edges[:-1],SM_NH[2], yerr=SM_NH_err[2], fmt='o-',label='Standard', markersize=5)
ax3.errorbar(bin_edges[:-1],FP_NH[2], yerr=FP_NH_err[2], fmt='o-',label='FP-NSSI', markersize=5)
ax3.errorbar(bin_edges[:-1],FV_NH[2], yerr=FV_NH_err[2], fmt='o-',label='FV-NSSI', markersize=5)

ax4 = fig2.add_subplot(234)
ax4.errorbar(bin_edges[:-1],SM_IH[0], yerr=SM_IH_err[0], fmt='o-',label='Standard', markersize=5)
ax4.errorbar(bin_edges[:-1],FP_IH[0], yerr=FP_IH_err[0], fmt='o-',label='FP-NSSI', markersize=5)
ax4.errorbar(bin_edges[:-1],FV_IH[0], yerr=FV_IH_err[0], fmt='o-',label='FV-NSSI', markersize=5)

ax5 = fig2.add_subplot(235)
ax5.errorbar(bin_edges[:-1],SM_IH[1], yerr=SM_IH_err[1], fmt='o-',label='Standard', markersize=5)
ax5.errorbar(bin_edges[:-1],FP_IH[1], yerr=FP_IH_err[1], fmt='o-',label='FP-NSSI', markersize=5)
ax5.errorbar(bin_edges[:-1],FV_IH[1], yerr=FV_IH_err[1], fmt='o-',label='FV-NSSI', markersize=5)

ax6 = fig2.add_subplot(236)
ax6.errorbar(bin_edges[:-1],SM_IH[2], yerr=SM_IH_err[2], fmt='o-',label='Standard', markersize=5)
ax6.errorbar(bin_edges[:-1],FP_IH[2], yerr=FP_IH_err[2], fmt='o-',label='FP-NSSI', markersize=5)
ax6.errorbar(bin_edges[:-1],FV_IH[2], yerr=FV_IH_err[2], fmt='o-',label='FV-NSSI', markersize=5)

ax1.set_ylabel('Events/1000ms')
ax1.set_xlabel('Time (ms)')
ax1.set_title('ES, NH', y=.9)
ax1.legend(loc='upper right');
ax1.set_xticks(xtick)
for tick in ax1.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax1.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax1.set_xlim(0,21000)

ax2.set_ylabel('Events/1000ms')
ax2.set_xlabel('Time (ms)')
ax2.set_title('IBD, NH', y=.9)
ax2.set_xticks(xtick)
for tick in ax2.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax2.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax2.set_xlim(0,21000)

ax3.set_ylabel('Events/1000ms')
ax3.set_xlabel('Time (ms)')
ax3.set_title('OSS, NH', y=.9)
ax3.set_xticks(xtick)
for tick in ax3.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax3.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax3.set_xlim(0,21000)

ax4.set_ylabel('Events/1000ms')
ax4.set_xlabel('Time (ms)')
ax4.set_title('ES, IH', y=.9)
ax4.set_xticks(xtick)
for tick in ax4.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax4.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax4.set_xlim(0,21000)

ax5.set_ylabel('Events/1000ms')
ax5.set_xlabel('Time (ms)')
ax5.set_title('IBD, IH', y=.9)
ax5.set_xticks(xtick)
for tick in ax5.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax5.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax5.set_xlim(0,21000)

ax6.set_ylabel('Events/1000ms')
ax6.set_xlabel('Time (ms)')
ax6.set_title('OSS, IH', y=.9)
ax6.set_xticks(xtick)
for tick in ax6.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax6.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax6.set_xlim(0,21000)

fig1.show()
fig2.show()
raw_input()

#plt.show()
