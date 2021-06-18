import numpy as np
import csv
import os
import matplotlib.pyplot as plt


NH = 'NH'
N_directories = ['osc000','osc010','osc001']
IH = 'IH'
I_directories = ['osc100','osc110','osc101']
all_times_sm_NH = np.array([])
all_times_g1_NH = np.array([])
all_times_g3_NH = np.array([])
all_times_sm_IH = np.array([])
all_times_g1_IH = np.array([])
all_times_g3_IH = np.array([])

def ratio_err(a,b,siga,sigb):
  return np.sqrt(float(siga*siga)/float(b*b) + float(sigb*sigb)*(float(a*a)/float(b**4)))

for num,d in enumerate(N_directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/nu_events/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith("1.txt"):# and ('es' in filepath):
        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          times = np.array(list(zip(*reader))[0], dtype=float)

        if (num == 0):
          all_times_sm_NH = np.concatenate([all_times_sm_NH,times], axis=None)
        if (num == 1):
          all_times_g1_NH = np.concatenate([all_times_g1_NH,times], axis=None)
        if(num == 2):
          all_times_g3_NH = np.concatenate([all_times_g3_NH,times], axis=None)

for num,d in enumerate(I_directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/nu_events/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith("1.txt"):# and ('es' in filepath):
        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          times = np.array(list(zip(*reader))[0], dtype=float)

        if (num == 0):
          all_times_sm_IH = np.concatenate([all_times_sm_IH,times], axis=None)
        if (num == 1):
          all_times_g1_IH = np.concatenate([all_times_g1_IH,times], axis=None)
        if(num == 2):
          all_times_g3_IH = np.concatenate([all_times_g3_IH,times], axis=None)

bins=range(1000,22000,1000)

all_times_sm_NH,bin_edges = np.histogram(all_times_sm_NH,bins=bins)
all_times_g1_NH,bin_edges = np.histogram(all_times_g1_NH,bins=bins)
all_times_g3_NH,bin_edges = np.histogram(all_times_g3_NH,bins=bins)

all_times_sm_IH,bin_edges = np.histogram(all_times_sm_IH,bins=bins)
all_times_g1_IH,bin_edges = np.histogram(all_times_g1_IH,bins=bins)
all_times_g3_IH,bin_edges = np.histogram(all_times_g3_IH,bins=bins)

error_sm_NH = np.sqrt(all_times_sm_NH)
error_g1_NH = np.sqrt(all_times_g1_NH)
error_g3_NH = np.sqrt(all_times_g3_NH)

error_sm_IH = np.sqrt(all_times_sm_IH)
error_g1_IH = np.sqrt(all_times_g1_IH)
error_g3_IH = np.sqrt(all_times_g3_IH)

SM_num_NH, SM_denom_NH, g1_num_NH, g1_denom_NH, g3_num_NH, g3_denom_NH = (0,)*6
SM_num_IH, SM_denom_IH, g1_num_IH, g1_denom_IH, g3_num_IH, g3_denom_IH = (0,)*6

for i in range(11,20):
  SM_num_NH += all_times_sm_NH[i]
  g1_num_NH += all_times_g1_NH[i]
  g3_num_NH += all_times_g3_NH[i]
  SM_num_IH += all_times_sm_IH[i]
  g1_num_IH += all_times_g1_IH[i]
  g3_num_IH += all_times_g3_IH[i]
for i in range(0,11):
  SM_denom_NH += all_times_sm_NH[i]
  g1_denom_NH += all_times_g1_NH[i]
  g3_denom_NH += all_times_g3_NH[i]
  SM_denom_IH += all_times_sm_IH[i]
  g1_denom_IH += all_times_g1_IH[i]
  g3_denom_IH += all_times_g3_IH[i]

fig = plt.figure()

ax2 = fig.add_subplot(121)
ax2.errorbar(bin_edges[:-1],all_times_sm_IH, yerr=error_sm_IH, fmt='o-', label='Standard', markersize=5)
ax2.errorbar(bin_edges[:-1],all_times_g3_IH, yerr=error_g3_IH, fmt='o-', label='FP-NSSI', markersize=5)
ax2.errorbar(bin_edges[:-1],all_times_g1_IH, yerr=error_g1_IH, fmt='o-', label='FV-NSSI', markersize=5)

ax1 = fig.add_subplot(122)
ax1.errorbar(bin_edges[:-1],all_times_sm_NH, yerr=error_sm_NH, fmt='o-', label='Standard', markersize=5)
ax1.errorbar(bin_edges[:-1],all_times_g3_NH, yerr=error_g3_NH, fmt='o-', label='FP-NSSI', markersize=5)
ax1.errorbar(bin_edges[:-1],all_times_g1_NH, yerr=error_g1_NH, fmt='o-', label='FV-NSSI', markersize=5)

print('NH Standard: R = %f +/- %f' %(float(SM_num_NH)/float(SM_denom_NH),ratio_err(SM_num_NH,SM_denom_NH,np.sqrt(SM_num_NH),np.sqrt(SM_denom_NH))))
print('NH FP-NSSI: R = %f +/- %f' %(float(g3_num_NH)/float(g3_denom_NH),ratio_err(g3_num_NH,g3_denom_NH,np.sqrt(g3_num_NH),np.sqrt(g3_denom_NH))))
print('NH FV-NSSI: R = %f +/- %f' %(float(g1_num_NH)/float(g1_denom_NH),ratio_err(g1_num_NH,g1_denom_NH,np.sqrt(g1_num_NH),np.sqrt(g1_denom_NH))))

print('IH Standard: R = %f +/- %f' %(float(SM_num_IH)/float(SM_denom_IH),ratio_err(SM_num_IH,SM_denom_IH,np.sqrt(SM_num_IH),np.sqrt(SM_denom_IH))))
print('IH FP-NSSI: R = %f +/- %f' %(float(g3_num_IH)/float(g3_denom_IH),ratio_err(g3_num_IH,g3_denom_IH,np.sqrt(g3_num_IH),np.sqrt(g3_denom_IH))))
print('IH FV-NSSI: R = %f +/- %f' %(float(g1_num_IH)/float(g1_denom_IH),ratio_err(g1_num_IH,g1_denom_IH,np.sqrt(g1_num_IH),np.sqrt(g1_denom_IH))))

xtick = [1000,3500,6000,8500,11000,13500,16000,18000,20000]

ax1.set_ylabel('Events/1000ms')
ax1.set_xlabel('Time (ms)')
ax1.set_title('All Events, '+NH, y=.9)
ax1.set_xticks(xtick)
for tick in ax1.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax1.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax1.set_xlim(0,21000)
ax1.legend(loc='upper right');

ax2.set_ylabel('Events/1000ms')
ax2.set_xlabel('Time (ms)')
ax2.set_title('All Events, '+IH, y=.9)
ax2.set_xticks(xtick)
for tick in ax2.xaxis.get_major_ticks():
  tick.label.set_fontsize(8)

for tick in ax2.yaxis.get_major_ticks():
  tick.label.set_fontsize(8)
ax2.set_xlim(0,21000)
ax2.legend(loc='upper right');

plt.show()
