import numpy as np
import csv
import os
import matplotlib.pyplot as plt



directories = ['osc100','osc110','osc101']
all_times_sm = np.array([])
all_times_g1 = np.array([])
all_times_g3 = np.array([])

for num,d in enumerate(directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/events_target/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith(".txt"):
        #print(filepath)

        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          first_col = list(zip(*reader))[0]

        times = np.array(first_col, dtype=float)
        #print(times)
        if (num == 0):
          all_times_sm = np.concatenate([all_times_sm,times], axis=None)
        if (num == 1):
          all_times_g1 = np.concatenate([all_times_sm,times], axis=None)
        if(num == 2):
          all_times_g3 = np.concatenate([all_times_sm,times], axis=None)

#print(times)
#print(all_times)

bins=range(1000,21000,1000)

all_times_sm,bin_edges = np.histogram(all_times_sm,bins=bins)
all_times_g1,bin_edges = np.histogram(all_times_g1,bins=bins)
all_times_g3,bin_edges = np.histogram(all_times_g3,bins=bins)

#print(len(all_times_hist))
#print(len(bin_edges[:-1]))

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(bin_edges[:-1],all_times_sm,c = 'b',label='SM')
ax1.scatter(bin_edges[:-1],all_times_g1,c = 'r',label='FV-NSSI')
ax1.scatter(bin_edges[:-1],all_times_g3,c = 'g',label='FP-NSSI')
#print(all_times_hist)
#print(bin_edges)
#plt.hist(times,bins=bins)
plt.ylabel('Number of events/1s')
plt.xlabel('time (ms)')
plt.xlim(800,20000)
plt.legend(loc='upper right');
plt.show()
